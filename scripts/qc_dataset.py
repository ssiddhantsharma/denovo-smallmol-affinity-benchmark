"""Deterministic integrity checks on the reference set. No network, no model, no GPU.

Runs in CI on every PR that touches reference/ or predictions/. It cannot invent
anything -- it only re-derives what the files already say and complains when they
disagree with each other. Every check is one a bad ingest has plausibly caused.

The two reference files are a pair: system_reference.csv says what a row IS
(protein, ligand, provenance) and experimental_reference_ground_truth.csv says
what was MEASURED. A row in one and not the other is the failure mode that
silently shrinks the benchmark, so that is checked first.

  python scripts/qc_dataset.py
  python scripts/qc_dataset.py --min-rows 60   # also fail if the set shrank
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SYSTEM = ROOT / "reference" / "system_reference.csv"
TRUTH = ROOT / "reference" / "experimental_reference_ground_truth.csv"
PREDICTIONS = ROOT / "predictions"

AA = set("ACDEFGHIKLMNPQRSTVWY")
SYSTEM_COLS = ["id", "protein", "ligand_name", "sequence", "smiles", "source", "source_url"]
TRUTH_COLS = ["id", "is_binder", "experimental_pKD"]
# pKd 2 is 10 mM, pKd 13 is 0.1 pM: outside this a unit slip is likelier than a measurement
PKD_RANGE = (2.0, 13.0)

checks: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    checks.append((name, bool(ok), detail))


def _rows(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh))


def main(min_rows: int) -> int:
    system, truth = _rows(SYSTEM), _rows(TRUTH)
    print(f"{SYSTEM.name}: {len(system)} rows | {TRUTH.name}: {len(truth)} rows")

    check("system_reference has the expected columns",
          list(system[0]) == SYSTEM_COLS, str(list(system[0])))
    check("ground_truth has the expected columns",
          list(truth[0]) == TRUTH_COLS, str(list(truth[0])))

    sids = [r["id"] for r in system]
    tids = [r["id"] for r in truth]
    check("ids are unique in system_reference", len(sids) == len(set(sids)),
          f"{len(sids) - len(set(sids))} duplicates")
    check("ids are unique in ground_truth", len(tids) == len(set(tids)),
          f"{len(tids) - len(set(tids))} duplicates")
    only_s, only_t = set(sids) - set(tids), set(tids) - set(sids)
    check("the two reference files cover the same ids", not only_s and not only_t,
          f"system-only {sorted(only_s)[:3]}, truth-only {sorted(only_t)[:3]}")

    check("row count did not shrink", len(system) >= min_rows,
          f"{len(system)} < {min_rows}")

    bad_seq = [r["id"] for r in system if not r["sequence"]
               or set(r["sequence"].upper()) - AA]
    check("sequences use the 20-letter alphabet", not bad_seq,
          f"{len(bad_seq)} rows, e.g. {bad_seq[:3]}")

    missing = [c for c in ("smiles", "source", "source_url", "ligand_name")
               for r in system if not r[c]]
    check("provenance and ligand fields are populated", not missing,
          f"{len(missing)} empty cells")

    # A pKd on a declared non-binder is contradictory; a binder without one cannot
    # enter the within-ligand correlation and would silently shrink it.
    nb_with_pkd, b_without = [], []
    for r in truth:
        has = (r["experimental_pKD"] or "").strip() not in ("", "nan")
        if r["is_binder"] == "1" and not has:
            b_without.append(r["id"])
        if r["is_binder"] == "0" and has:
            nb_with_pkd.append(r["id"])
    check("every binder has a pKd", not b_without, f"{b_without[:3]}")
    check("no non-binder carries a pKd", not nb_with_pkd, f"{nb_with_pkd[:3]}")
    check("is_binder is 0 or 1", all(r["is_binder"] in ("0", "1") for r in truth),
          str(sorted({r["is_binder"] for r in truth})))

    pkds = [float(r["experimental_pKD"]) for r in truth
            if (r["experimental_pKD"] or "").strip() not in ("", "nan")]
    out_of_range = [p for p in pkds if not (PKD_RANGE[0] <= p <= PKD_RANGE[1])]
    check(f"pKd within {PKD_RANGE[0]}-{PKD_RANGE[1]}", not out_of_range,
          f"{len(out_of_range)} outside, e.g. {out_of_range[:3]}")
    check("pKd values are finite", all(math.isfinite(p) for p in pkds))

    # The within-ligand statistic needs at least two binders sharing a ligand;
    # a ligand with one binder contributes nothing and should be visible, not silent.
    per_ligand: dict[str, int] = {}
    binder_ids = {r["id"] for r in truth if r["is_binder"] == "1"}
    for r in system:
        if r["id"] in binder_ids:
            per_ligand[r["ligand_name"]] = per_ligand.get(r["ligand_name"], 0) + 1
    usable = {k: v for k, v in per_ligand.items() if v >= 2}
    print(f"  ligands with >=2 binders (these carry the within-ligand result): "
          f"{len(usable)} of {len(per_ligand)}")
    check("at least 3 ligands support the within-ligand statistic", len(usable) >= 3,
          f"{len(usable)} ligands with 2+ binders")

    # Every predictions file must cover exactly the reference ids, or score.py is
    # quietly averaging over a different set for that method than for the others.
    pred_files = sorted(PREDICTIONS.glob("*_predictions.csv"))
    check("prediction files exist", bool(pred_files), f"{len(pred_files)} files")
    mismatched = []
    for p in pred_files:
        pid = {r["id"] for r in _rows(p)}
        if pid != set(sids):
            mismatched.append(f"{p.name} ({len(set(sids) - pid)} missing, "
                              f"{len(pid - set(sids))} extra)")
    check("every predictions file covers exactly the reference ids", not mismatched,
          "; ".join(mismatched[:3]))

    # rdkit is an optional dependency locally but installed in CI, so a missing
    # import downgrades this check rather than failing the run
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        bad = [r["id"] for r in system if Chem.MolFromSmiles(r["smiles"]) is None]
        check("every SMILES parses", not bad, f"{len(bad)} rows, e.g. {bad[:3]}")
    except ImportError:
        print("  (rdkit not installed: SMILES parsing not checked)")

    width = max(len(n) for n, _, _ in checks)
    print()
    for name, ok, detail in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:{width}}  {detail if not ok else ''}")
    failed = sum(1 for _, ok, _ in checks if not ok)
    print(f"\n{len(checks) - failed}/{len(checks)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-rows", type=int, default=0)
    raise SystemExit(main(ap.parse_args().min_rows))
