"""Derive interface-PAE metrics from saved fold outputs -> predictions/*.csv. Run after
run_predictors boltz+protenix (needs the _work/ outputs). Lower PAE = tighter interface.

Boltz: interface block of the full PAE matrix (protein residues 0..L, ligand atoms L..N; L = len
of the protein sequence). Protenix: chain_pair_pae_{min,mean}[protein][ligand] from the summary
JSON, averaged over samples. Needs numpy.
"""

import csv
import glob
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "_work"


def refs():
    return list(csv.DictReader((ROOT / "reference/system_reference.csv").read_text().splitlines()))


def write(method, vals):
    with open(ROOT / f"predictions/{method}_predictions.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["method", "id", "predicted_affinity"])
        for r in refs():
            v = vals.get(r["id"])
            w.writerow([method, r["id"], "" if v is None else v])
    print(f"wrote {method}_predictions.csv  (n={len(vals)})")


def boltz():
    pred = WORK / "boltz_out/boltz_results_boltz_yaml/predictions"
    pmin, pmean = {}, {}
    for i, r in enumerate(refs()):
        hits = glob.glob(str(pred / f"c{i:02d}" / "pae*.npz"))
        if not hits:
            continue
        pae = np.load(hits[0])["pae"]; lp = len(r["sequence"])
        block = np.concatenate([pae[:lp, lp:].ravel(), pae[lp:, :lp].ravel()])
        pmin[r["id"]] = float(block.min()); pmean[r["id"]] = float(block.mean())
    write("boltz-2-pae-min", pmin); write("boltz-2-pae-mean", pmean)


def protenix():
    pmin, pmean = {}, {}
    for i, r in enumerate(refs()):
        mn, me = [], []
        for f in sorted(glob.glob(str(WORK / f"px_out/d{i}/seed_101/predictions/"
                                        f"d{i}_summary_confidence_sample_*.json"))):
            d = json.loads(Path(f).read_text())
            mn.append(d["chain_pair_pae_min"][0][1]); me.append(d["chain_pair_pae_mean"][0][1])
        if mn:
            pmin[r["id"]] = sum(mn) / len(mn); pmean[r["id"]] = sum(me) / len(me)
    write("protenix-pae-min", pmin); write("protenix-pae-mean", pmean)


if __name__ == "__main__":
    boltz()
    protenix()
