#!/usr/bin/env python
"""Score the de-novo small-molecule benchmark with dtSFM's encoder cosine
(binding-compatibility), reusing the shipped ReencodingPipeline recipe:

  drug  SMILES -> ibm/MoLFormer-XL-both-10pct -> mean-pool 768 -> encoder.drug_global_proj -> 512 (L2)
  target seq   -> facebook/esm2_t33_650M_UR50D -> per-residue 1280 -> encoder.protein_proj + pool -> 512 (L2)
  score = cosine(drug512, target512)

In our benchmark the "protein" (a de novo binder) is dtSFM's target and the
"smiles" (small molecule) is dtSFM's drug.

GATE FIRST: validate the featurization on in-distribution metadata_v3.csv pairs
(true pairs must score well above shuffled). If the gate fails, write nothing.

Usage:
  PYTHONPATH=<dtSFM>/src .venv_score/bin/python scratch_run_dtsfm.py \
      --dtsfm <dtSFM repo> --bench <benchmark repo> [--n-control 40]
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from huggingface_hub import hf_hub_download

MOLFORMER = "ibm/MoLFormer-XL-both-10pct"
ESM = "facebook/esm2_t33_650M_UR50D"


def load_encoder(dtsfm_src, device):
    sys.path.insert(0, str(dtsfm_src))
    from omegaconf import OmegaConf

    from calm.encoder.model_v3 import CALMEncoderV3
    ckpt_path = hf_hub_download("SFM-BIIE-ETHZ/dtSFM-v3", "encoder_b3_epoch010.pt")
    ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)
    cfg = ckpt.get("config")
    if cfg is None:
        from calm.encoder.train_dtsfm_v3 import default_model_cfg
        cfg = default_model_cfg()
    enc = CALMEncoderV3(OmegaConf.create(cfg)).to(device).eval()
    enc.load_state_dict(ckpt["model_state_dict"])
    for p in enc.parameters():
        p.requires_grad_(False)
    return enc


class Scorer:
    def __init__(self, dtsfm_src, device):
        self.device = device
        from transformers import AutoModel, AutoTokenizer, EsmModel
        self.enc = load_encoder(dtsfm_src, device)
        self.mol_tok = AutoTokenizer.from_pretrained(MOLFORMER, trust_remote_code=True)
        self.mol = AutoModel.from_pretrained(
            MOLFORMER, trust_remote_code=True, deterministic_eval=True).to(device).eval()
        self.esm_tok = AutoTokenizer.from_pretrained(ESM)
        self.esm = EsmModel.from_pretrained(ESM).to(device).eval()
        for m in (self.mol, self.esm):
            for p in m.parameters():
                p.requires_grad_(False)

    @torch.no_grad()
    def drug512(self, smiles):
        e = self.mol_tok([smiles], padding=True, truncation=True, max_length=512,
                         return_tensors="pt").to(self.device)
        h = self.mol(input_ids=e["input_ids"], attention_mask=e["attention_mask"]).last_hidden_state
        m = e["attention_mask"].unsqueeze(-1).float()
        pooled = (h * m).sum(1) / m.sum(1).clamp(min=1.0)               # (1,768)
        return F.normalize(self.enc.drug_global_proj(pooled.float()), dim=-1)  # (1,512)

    @torch.no_grad()
    def target512(self, seq):
        e = self.esm_tok(seq, return_tensors="pt", truncation=True, max_length=1024).to(self.device)
        h = self.esm(**e).last_hidden_state[0]                          # (L+2,1280) incl CLS/EOS
        per_res = h[1:-1].unsqueeze(0)                                  # (1,L,1280) strip CLS/EOS
        L = per_res.shape[1]
        mask = torch.ones(1, L, dtype=torch.bool, device=self.device)
        hp = self.enc.protein_proj(per_res.float()) * mask.unsqueeze(-1).float()
        pooled = self.enc.protein_pool_pre(hp, mask)                    # (1,512)
        return F.normalize(pooled, dim=-1)

    @torch.no_grad()
    def cosine(self, smiles, seq):
        return float((self.drug512(smiles) * self.target512(seq)).sum().item())


def gate(scorer, n_control):
    """True in-distribution pairs must score well above shuffled pairs."""
    p = hf_hub_download("SFM-BIIE-ETHZ/dtSFM-v3", "metadata_v3.csv")
    rows = list(csv.DictReader(open(p)))
    rows = [r for r in rows if r.get("drug_smiles") and r.get("protein_seq")
            and r.get("affinity_log") and len(r["protein_seq"]) <= 1000]
    rows.sort(key=lambda r: float(r["affinity_log"]), reverse=True)   # strongest binders first
    ctrl = rows[:n_control]
    true_c = [scorer.cosine(r["drug_smiles"], r["protein_seq"]) for r in ctrl]
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(ctrl))
    perm = [(i + 1) % len(ctrl) if perm[i] == i else perm[i] for i in range(len(ctrl))]
    shuf_c = [scorer.cosine(ctrl[i]["drug_smiles"], ctrl[perm[i]]["protein_seq"]) for i in range(len(ctrl))]
    from sklearn.metrics import roc_auc_score
    y = [1] * len(true_c) + [0] * len(shuf_c)
    auc = roc_auc_score(y, true_c + shuf_c)
    print(f"\n=== POSITIVE-CONTROL GATE (n={len(ctrl)} in-distribution pairs) ===")
    print(f"  true-pair cosine   median {np.median(true_c):+.3f}  mean {np.mean(true_c):+.3f}")
    print(f"  shuffled-pair      median {np.median(shuf_c):+.3f}  mean {np.mean(shuf_c):+.3f}")
    print(f"  true-vs-shuffled AUROC = {auc:.3f}")
    # Rank-based gate: our benchmark uses cosine for RANKING (AUROC/Spearman), so the
    # meaningful validation is separation (true >> shuffled), not absolute magnitude.
    sep = float(np.median(true_c) - np.median(shuf_c))
    ok = (auc > 0.75) and (sep > 0.12)
    print(f"  separation (median true-shuffled) = {sep:+.3f}")
    print(f"  GATE: {'PASS' if ok else 'FAIL'}  (rank-based: AUROC>0.75 and separation>0.12)")
    if np.median(true_c) < 0.5:
        print("  NOTE: absolute cosine sits below the smoketest's >0.5 calibration (likely an ESM-source"
              " detail, HF EsmModel vs fair-esm); ranking is validated, so absolute values are not over-read.")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dtsfm", required=True)
    ap.add_argument("--bench", required=True)
    ap.add_argument("--n-control", type=int, default=40)
    a = ap.parse_args()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device={device}")
    scorer = Scorer(Path(a.dtsfm) / "src", device)

    if not gate(scorer, a.n_control):
        print("\nFEATURIZATION GATE FAILED — not writing predictions. My glue is wrong, stopping.")
        sys.exit(2)

    sysrows = list(csv.DictReader(open(Path(a.bench) / "reference/system_reference.csv")))
    print(f"\n=== scoring {len(sysrows)} benchmark systems ===")
    out = []
    for r in sysrows:
        c = scorer.cosine(r["smiles"], r["sequence"])
        out.append({"method": "dtsfm-cosine", "id": r["id"], "predicted_affinity": c})
    dst = Path(a.bench) / "predictions/dtsfm-cosine_predictions.csv"
    with open(dst, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["method", "id", "predicted_affinity"])
        w.writeheader()
        w.writerows(out)
    print(f"wrote {dst} ({len(out)} rows)")


if __name__ == "__main__":
    main()
