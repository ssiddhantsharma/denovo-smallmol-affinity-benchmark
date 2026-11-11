"""Deterministic ingest of one paper's extracted rows into the reference set.

This is the trustworthy half of the pipeline. An agent produces a flat CSV of
what it read; THIS script decides what enters the benchmark. No model writes
reference/ -- so a hallucinated row cannot appear without passing checks that
only look at values already in the extract.

Expected extract columns (header-flexible, missing -> empty):
  protein, ligand_name, sequence, smiles, kd_nM, is_binder, source, source_url, notes

What it enforces, and why each one has bitten a dataset like this:
  * an id is DERIVED, never taken from the extract, so two papers cannot collide
    on a name the model chose;
  * a row needs a parseable SMILES and a 20-letter sequence, or the predictors
    cannot score it and it would sit in the set as a permanent gap;
  * a binder needs a KD and a non-binder must not have one;
  * dedup is on (sequence, canonical SMILES), so the same pair re-extracted from
    a second paper does not enter twice under a different name;
  * NEW ROWS HAVE NO PREDICTIONS. They are written to reference/ and listed, and
    score.py will refuse them until the predictors are run, which needs a GPU and
    is not part of this step.

  python scripts/append_paper.py --extract extract.csv --dry-run
  python scripts/append_paper.py --extract extract.csv
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SYSTEM = ROOT / "reference" / "system_reference.csv"
TRUTH = ROOT / "reference" / "experimental_reference_ground_truth.csv"
BACKUPS = ROOT / "reference" / "backups"

AA = set("ACDEFGHIKLMNPQRSTVWY")
SYSTEM_COLS = ["id", "protein", "ligand_name", "sequence", "smiles", "source", "source_url"]
TRUTH_COLS = ["id", "is_binder", "experimental_pKD"]
PKD_RANGE = (2.0, 13.0)


def canonical_smiles(smiles: str) -> str | None:
    try:
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
    except ImportError as exc:                     # refuse rather than dedup on raw text
        raise SystemExit("rdkit is required to ingest: pip install rdkit") from exc
    mol = Chem.MolFromSmiles((smiles or "").strip())
    return Chem.MolToSmiles(mol) if mol is not None else None


def slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", (text or "").lower()).strip("-")[:40]


def read_rows(path: Path) -> list[dict]:
    with path.open() as fh:
        return [{(k or "").strip(): (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(fh)]


def main(extract: Path, dry_run: bool) -> None:
    rows = read_rows(extract)
    system = read_rows(SYSTEM)

    existing_pairs = {}
    for r in system:
        cs = canonical_smiles(r["smiles"])
        existing_pairs[(r["sequence"], cs)] = r["id"]
    taken = {r["id"] for r in system}

    fresh_system, fresh_truth, rejected = [], [], []
    for i, r in enumerate(rows, 1):
        why = None
        seq = (r.get("sequence") or "").upper()
        cs = canonical_smiles(r.get("smiles", ""))
        is_binder = (r.get("is_binder") or "").strip().lower() in ("1", "true", "yes", "y")
        kd = (r.get("kd_nM") or "").strip()

        if not seq or set(seq) - AA:
            why = "sequence missing or not the 20-letter alphabet"
        elif cs is None:
            why = f"SMILES does not parse: {r.get('smiles', '')[:40]!r}"
        elif not r.get("source") or not r.get("source_url"):
            why = "no source / source_url"
        elif is_binder and not kd:
            why = "binder with no KD"
        elif not is_binder and kd:
            why = "non-binder carrying a KD"
        elif (seq, cs) in existing_pairs:
            why = f"already present as {existing_pairs[(seq, cs)]}"

        pkd = ""
        if why is None and is_binder:
            try:
                pkd_val = -__import__("math").log10(float(kd) * 1e-9)
            except (ValueError, OverflowError):
                why = f"KD is not a positive number: {kd!r}"
            else:
                if not (PKD_RANGE[0] <= pkd_val <= PKD_RANGE[1]):
                    why = f"pKd {pkd_val:.2f} outside {PKD_RANGE} - unit slip?"
                else:
                    pkd = repr(pkd_val)

        if why:
            rejected.append((i, r.get("protein", "?"), r.get("ligand_name", "?"), why))
            continue

        base = f"{slug(r.get('protein', 'design'))}__{slug(r.get('ligand_name', 'ligand'))}"
        new_id, n = base, 2
        while new_id in taken:
            new_id, n = f"{base}-{n}", n + 1
        taken.add(new_id)
        existing_pairs[(seq, cs)] = new_id

        fresh_system.append({"id": new_id, "protein": r.get("protein", ""),
                             "ligand_name": r.get("ligand_name", ""), "sequence": seq,
                             "smiles": r.get("smiles", "").strip(),
                             "source": r.get("source", ""), "source_url": r.get("source_url", "")})
        fresh_truth.append({"id": new_id, "is_binder": "1" if is_binder else "0",
                            "experimental_pKD": pkd})

    print(f"extract: {len(rows)} rows | new: {len(fresh_system)} | rejected: {len(rejected)}")
    for i, prot, lig, why in rejected:
        print(f"  row {i} [{prot} / {lig}]: {why}")
    for s in fresh_system:
        print(f"  + {s['id']}")

    if dry_run:
        print("[dry-run] nothing written.")
        return
    if not fresh_system:
        print("nothing to add.")
        return

    BACKUPS.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    for path in (SYSTEM, TRUTH):
        shutil.copy(path, BACKUPS / f"{path.stem}_{stamp}.csv")

    with SYSTEM.open("a", newline="") as fh:
        csv.DictWriter(fh, fieldnames=SYSTEM_COLS).writerows(fresh_system)
    with TRUTH.open("a", newline="") as fh:
        csv.DictWriter(fh, fieldnames=TRUTH_COLS).writerows(fresh_truth)

    print(f"backed up to {BACKUPS.relative_to(ROOT)}/, appended {len(fresh_system)} rows")
    print("\nNEXT: these rows have no predictions yet. Until the predictors are run "
          "they are in reference/ but not in any figure, and scripts/qc_dataset.py "
          "will fail the 'predictions cover exactly the reference ids' check -- which "
          "is the intended reminder, not a bug. See 'Adding a paper' in the README.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", required=True, type=Path)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    main(a.extract, a.dry_run)
