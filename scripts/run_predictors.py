"""Run each predictor over reference/system_reference.csv; write predictions/<method>_predictions.csv
(strict method,id,predicted_affinity). Analysis lives in score.py. Cofolders need their own GPU envs.

  boltz      boltz CLI, affinity, single-sequence (no MSA), recycling 3 -> boltz-2{,-pbind,-iptm}
  nesso      recursionpharma/nesso, seq+SMILES YAML (ESM)             -> nesso-1{,-pbind}
  protenix   Protenix-v2, single-sequence, seed 101                   -> protenix-{iptm,ligiptm,gpde,ranking}
  baselines  RDKit MW and cLogP (no GPU)                              -> molecular-weight, clogp

Env: CUDA_VISIBLE_DEVICES, BOLTZ, NESSO, PROTENIX_DIR, PROTENIX_PY.
Cofolders emit affinity_pred_value; OpenBind pK convention: pK = 6 - affinity_pred_value.
"""

import argparse
import csv
import json
import os
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "predictions"
WORK = ROOT / "_work"; WORK.mkdir(exist_ok=True)
GPU = os.environ.get("CUDA_VISIBLE_DEVICES", "0")
BOLTZ = os.environ.get("BOLTZ", "boltz")
NESSO = os.environ.get("NESSO", "nesso")


def refs():
    return list(csv.DictReader((ROOT / "reference/system_reference.csv").read_text().splitlines()))


def write_pred(method, values):
    """values: {id -> float or None}."""
    with open(PRED / f"{method}_predictions.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["method", "id", "predicted_affinity"])
        for r in refs():
            v = values.get(r["id"])
            w.writerow([method, r["id"], "" if v is None else v])
    print(f"WROTE {method}_predictions.csv", flush=True)


def yaml_for(seq, smiles, msa_empty):
    msa = "      msa: empty\n" if msa_empty else ""
    return (f"sequences:\n  - protein:\n      id: A\n      sequence: {seq}\n{msa}"
            f"  - ligand:\n      id: B\n      smiles: '{smiles}'\n"
            f"properties:\n  - affinity:\n      binder: B\n")


def run_boltz():
    yd = WORK / "boltz_yaml"; yd.mkdir(exist_ok=True)
    for i, r in enumerate(refs()):
        (yd / f"c{i:02d}.yaml").write_text(yaml_for(r["sequence"], r["smiles"], True))
    subprocess.run([BOLTZ, "predict", str(yd), "--out_dir", str(WORK / "boltz_out"),
                    "--recycling_steps", "3", "--diffusion_samples", "1", "--seed", "101",
                    "--override", "--no_kernels", "--write_full_pae", "--devices", "1"],
                   env={**os.environ, "CUDA_VISIBLE_DEVICES": GPU}, check=True)
    p = WORK / "boltz_out" / "boltz_results_boltz_yaml" / "predictions"
    pk, pb, it = {}, {}, {}
    for i, r in enumerate(refs()):
        aff = p / f"c{i:02d}" / f"affinity_c{i:02d}.json"
        cf = p / f"c{i:02d}" / f"confidence_c{i:02d}_model_0.json"
        if not aff.exists():
            print(f"MISSING boltz {r['id']}", flush=True); continue
        a = json.loads(aff.read_text()); c = json.loads(cf.read_text()) if cf.exists() else {}
        v = a.get("affinity_pred_value")
        pk[r["id"]] = 6 - v if v is not None else None
        pb[r["id"]] = a.get("affinity_probability_binary"); it[r["id"]] = c.get("iptm")
    write_pred("boltz-2", pk); write_pred("boltz-2-pbind", pb); write_pred("boltz-2-iptm", it)


def run_nesso():
    yd = WORK / "nesso_yaml"; yd.mkdir(exist_ok=True)
    for i, r in enumerate(refs()):
        (yd / f"c{i:02d}.yaml").write_text(yaml_for(r["sequence"], r["smiles"], False))
    subprocess.run([NESSO, "predict", str(yd), "--out_dir", str(WORK / "nesso_out"),
                    "--override", "--devices", "1"],
                   env={**os.environ, "CUDA_VISIBLE_DEVICES": GPU}, check=True)
    pk, pb = {}, {}
    for i, r in enumerate(refs()):
        aff = WORK / "nesso_out" / "predictions" / f"c{i:02d}" / "affinity.json"
        if not aff.exists():
            print(f"MISSING nesso {r['id']}", flush=True); continue
        a = json.loads(aff.read_text()); v = a.get("affinity_pred_value")
        pk[r["id"]] = 6 - v if v is not None else None
        pb[r["id"]] = a.get("affinity_probability_binary")
    write_pred("nesso-1", pk); write_pred("nesso-1-pbind", pb)


def run_protenix():
    pdir = os.environ["PROTENIX_DIR"]; py = os.environ.get("PROTENIX_PY", "python")
    J = WORK / "protenix.json"
    J.write_text(json.dumps([{"seq": r["sequence"], "ligand": r["smiles"], "id": r["id"]}
                             for r in refs()], indent=2))
    sc = str(ROOT / "scripts/protenix_score.py")
    subprocess.run([py, sc, "build", "--in", str(J), "--input-json", str(WORK / "px_in.json")], check=True)
    cmd = (f'cd "{pdir}" && CUDA_VISIBLE_DEVICES={GPU} .venv/bin/protenix pred '
           f'-i {WORK}/px_in.json -o {WORK}/px_out -s 101 -n protenix-v2 '
           f'--use_msa false --use_default_params true')
    subprocess.run(["bash", "-c", cmd], check=True)
    subprocess.run([py, sc, "parse", "--in", str(J), "--outdir", str(WORK / "px_out"), "--seed", "101"],
                   check=True)
    s = {r["id"]: r for r in json.loads(J.read_text())}
    for method, col in [("protenix-iptm", "protenix_iptm"), ("protenix-ligiptm", "protenix_lig_iptm"),
                        ("protenix-gpde", "protenix_gpde"), ("protenix-ranking", "protenix_ranking")]:
        write_pred(method, {i: s[i].get(col) for i in s})


def run_baselines():
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors
    for method, f in [("clogp", Crippen.MolLogP), ("molecular-weight", Descriptors.MolWt)]:
        vals = {}
        for r in refs():
            m = Chem.MolFromSmiles(r["smiles"]); vals[r["id"]] = f(m) if m else None
        write_pred(method, vals)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["boltz", "nesso", "protenix", "baselines"])
    {"boltz": run_boltz, "nesso": run_nesso, "protenix": run_protenix,
     "baselines": run_baselines}[ap.parse_args().mode]()


if __name__ == "__main__":
    main()
