"""Extra Protenix confidence metrics from the saved summaries (no GPU, no re-fold).

pTM and per-chain pLDDT are in the Anthropic de-novo metric set but were not emitted by the
original Protenix predictor. They are already present in the saved summary jsons, so this reads the
best-ranked sample per system (no re-folding) and writes them as predictions.

  protenix-ptm            global pTM
  protenix-plddt-binder   chain A (binder protein) pLDDT
  protenix-plddt-ligand   chain B (ligand) pLDDT

ipSAE and LIS (Anthropic's strongest metrics) are protein-protein interface scores built from the
full residue-residue PAE matrix, which Protenix did not save here; their protein-small-molecule
analog is lig-ipTM + interface-PAE, already in the benchmark. Higher is more confident.
"""

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "predictions"
PX_OUT = ROOT / "_work" / "px_out"
SEED = 101


def refs():
    return list(csv.DictReader((ROOT / "reference/system_reference.csv").read_text().splitlines()))


def best_summary(i):
    d = PX_OUT / f"d{i}" / f"seed_{SEED}" / "predictions"
    best, best_rank = None, -1.0
    for j in sorted(d.glob(f"d{i}_summary_confidence_sample_*.json")):
        try:
            s = json.loads(j.read_text())
        except (FileNotFoundError, ValueError):
            continue
        if (s.get("ranking_score", -1.0) or -1.0) > best_rank:
            best, best_rank = s, s.get("ranking_score", -1.0)
    return best


def main():
    ptm, pb, pl = {}, {}, {}
    for i, r in enumerate(refs()):
        s = best_summary(i)
        if not s:
            print(f"  MISSING summary for {r['id']} (d{i})", flush=True)
            continue
        ptm[r["id"]] = s.get("ptm")
        cp = s.get("chain_plddt") or []
        if len(cp) >= 2:
            pb[r["id"]], pl[r["id"]] = cp[0], cp[1]
    for method, vals in [("protenix-ptm", ptm), ("protenix-plddt-binder", pb),
                         ("protenix-plddt-ligand", pl)]:
        with open(PRED / f"{method}_predictions.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["method", "id", "predicted_affinity"])
            for r in refs():
                v = vals.get(r["id"])
                w.writerow([method, r["id"], "" if v is None else v])
        print(f"WROTE {method}_predictions.csv ({len(vals)}/61)", flush=True)


if __name__ == "__main__":
    main()
