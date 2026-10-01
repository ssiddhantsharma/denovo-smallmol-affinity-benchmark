"""Pose-derived interface metrics from the Protenix co-folds (no GPU, no network).

These read the interface *geometry* rather than ligand composition, so unlike cLogP they are not
constant within a ligand and can in principle rank the de-novo proteins that bind one molecule:

  interface-ncontacts   protein-ligand heavy-atom pairs within 4.5 A (interface size)
  interface-ligburial   fraction of ligand heavy atoms with a protein heavy atom within 4.5 A (burial)

Reads the best-ranked Protenix sample from _work/px_out (produced by run_predictors.py protenix) and
writes predictions/interface-{ncontacts,ligburial}_predictions.csv. Higher is more buried/contacted.
"""

import csv
import json
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "predictions"
PX_OUT = ROOT / "_work" / "px_out"
SEED = 101
CUT = 4.5


def refs():
    return list(csv.DictReader((ROOT / "reference/system_reference.csv").read_text().splitlines()))


def best_pose(i):
    d = PX_OUT / f"d{i}" / f"seed_{SEED}" / "predictions"
    best, best_rank = None, -1.0
    for cif in sorted(d.glob(f"d{i}_sample_*.cif")):
        n = cif.stem.split("_sample_")[1]
        try:
            rk = json.loads((d / f"d{i}_summary_confidence_sample_{n}.json").read_text()).get("ranking_score", -1.0)
        except (FileNotFoundError, ValueError):
            rk = -1.0
        if rk > best_rank:
            best, best_rank = cif, rk
    return best


def heavy_coords(residues):
    return np.array([[a.pos.x, a.pos.y, a.pos.z]
                     for r in residues for a in r if a.element.name != "H"])


def interface(cif):
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    model = st[0]
    prot, lig = [], []
    for ch in model:
        poly = ch.get_polymer()
        (prot if poly and len(poly) > 10 else lig).extend(ch)
    P, L = heavy_coords(prot), heavy_coords(lig)
    if not len(P) or not len(L):
        return None
    dmin = np.sqrt(((L[:, None] - P[None]) ** 2).sum(-1))  # (lig, prot)
    ncontacts = int((dmin < CUT).sum())
    ligburial = float((dmin.min(1) < CUT).mean())
    return ncontacts, ligburial


def main():
    ncontacts, ligburial = {}, {}
    for i, r in enumerate(refs()):
        cif = best_pose(i)
        res = interface(cif) if cif and cif.exists() else None
        if res is None:
            print(f"  MISSING pose for {r['id']} (d{i})", flush=True)
            continue
        ncontacts[r["id"]], ligburial[r["id"]] = res
    for method, vals in [("interface-ncontacts", ncontacts), ("interface-ligburial", ligburial)]:
        with open(PRED / f"{method}_predictions.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["method", "id", "predicted_affinity"])
            for r in refs():
                v = vals.get(r["id"])
                w.writerow([method, r["id"], "" if v is None else v])
        print(f"WROTE {method}_predictions.csv ({len(vals)}/61)", flush=True)


if __name__ == "__main__":
    main()
