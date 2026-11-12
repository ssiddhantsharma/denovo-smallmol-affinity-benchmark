"""Bring the benchmark up to date: score whatever has no predictions yet.

The one command to run on the GPU box after `append_paper.py` has added rows.
It finds reference rows that no predictions file covers, co-folds them, scores
them, and leaves the rest to CI -- pushing `predictions/` triggers the figures
workflow, which re-scores and redraws.

Two ways to get predictions, and only one of them needs this box:

  scripts/run_rowan_cofold.py   hosted end to end. Rowan folds the complex from
                                sequence + SMILES and returns Boltz-2's affinity
                                with it, so no GPU and no local pose. Prefer this.
  --cofold + --rowan (below)    the local path: co-fold here, then submit the
                                pose to Rowan's pose-based rescorers (SQM,
                                GNINA, AEV-PLIG), which take a pose as input and
                                cannot start from a sequence.

So the GPU is needed for the pose-based rescorers, not for scoring as such. An
earlier version of this file claimed the GPU step was unavoidable; that was
wrong, and `rowan.submit_protein_cofolding_workflow` is why.

What stays off a public repo's Actions is the API KEY, not the compute: Rowan
charges credits, and `--max-credits` caps a workflow rather than a month.

It spends money, so it is dry-run by default and every spending step is capped:

  python scripts/update.py                     # what WOULD run, no GPU, no credits
  python scripts/update.py --cofold            # the GPU half only
  python scripts/update.py --rowan --max-credits 200    # the paid half, capped
  python scripts/update.py --all --max-credits 200

Credits are per workflow in Rowan's API, not per run, so --max-credits alone
cannot bound a session. `--max-systems` does, and is the knob to reach for when
a paper lands with forty pairs.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

import methods

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "system_reference.csv"
PRED = ROOT / "predictions"
WORK = ROOT / "_work"


def reference_ids() -> list[str]:
    with REF.open() as fh:
        return [r["id"] for r in csv.DictReader(fh)]


def covered_ids() -> dict[str, set[str]]:
    """method -> ids it has a prediction for."""
    out = {}
    for path in sorted(PRED.glob("*_predictions.csv")):
        with path.open() as fh:
            out[path.stem.replace("_predictions", "")] = {r["id"] for r in csv.DictReader(fh)}
    return out


def absent(covered: dict[str, set[str]]) -> list[str]:
    """Registered methods with no predictions file at all.

    These are invisible to `missing()`, which can only compare files that exist.
    That is how this script once reported "nothing to update" while three
    hosted co-folding columns had never been produced.
    """
    return [m for m in methods.expected() if m not in covered]


def missing(ids: list[str], covered: dict[str, set[str]]) -> dict[str, list[str]]:
    """method -> reference ids it is missing. A method missing EVERY id is a
    method that was never run; one missing a few is the case this script exists
    for, and the two want different responses, so they are reported apart."""
    return {m: [i for i in ids if i not in have] for m, have in covered.items()}


def _run(cmd: list[str], env_note: str = "") -> int:
    print(f"\n$ {' '.join(cmd)}{('   # ' + env_note) if env_note else ''}")
    return subprocess.call(cmd)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cofold", action="store_true", help="run the GPU co-folding step")
    ap.add_argument("--rowan", action="store_true", help="submit to Rowan (SPENDS CREDITS)")
    ap.add_argument("--all", action="store_true", help="both halves")
    ap.add_argument("--max-credits", type=int, default=None,
                    help="per-workflow credit cap passed to run_rowan.py")
    ap.add_argument("--max-systems", type=int, default=None,
                    help="refuse to start if more than this many systems need scoring")
    ap.add_argument("--poses", default="protenix", choices=("protenix", "boltz"))
    ap.add_argument("--count", action="store_true",
                    help="print only how many systems need scoring, for cron. Parsing "
                         "the prose output is how the weekly job silently stopped "
                         "submitting when a column was renamed.")
    a = ap.parse_args()
    do_cofold, do_rowan = a.cofold or a.all, a.rowan or a.all

    ids = reference_ids()
    covered = covered_ids()
    gaps = missing(ids, covered)
    never_produced = absent(covered)
    never_run = sorted(m for m, g in gaps.items() if len(g) == len(ids))
    partial = {m: g for m, g in gaps.items() if g and len(g) < len(ids)}
    need = sorted({i for g in partial.values() for i in g})

    if a.count:
        print(len(need))
        return 0
    print(f"reference: {len(ids)} systems | predictions files: {len(covered)}")
    if never_produced:
        print(f"\nregistered methods with NO predictions file ({len(never_produced)}):")
        for m in never_produced:
            flag = "optional" if m in methods.OPTIONAL else "MISSING"
            print(f"  [{flag}] {m:28s} <- {methods.PRODUCERS[m]}")
    if never_run:
        print(f"methods with no predictions at all ({len(never_run)}): {', '.join(never_run)}")
        print("  these are a full re-run, not an update; this script will not start one")
    if not partial:
        print("\nevery predictions file that exists covers every reference row.")
        if not never_produced:
            print("nothing to update.")
        return 0
    print(f"\nsystems needing scores: {len(need)}")
    for m, g in sorted(partial.items()):
        print(f"  {m:28s} missing {len(g)}: {', '.join(g[:4])}{' ...' if len(g) > 4 else ''}")

    if a.max_systems is not None and len(need) > a.max_systems:
        print(f"\nrefusing: {len(need)} systems exceeds --max-systems {a.max_systems}. "
              "Raise it deliberately, or ingest fewer rows at a time.")
        return 1

    if not (do_cofold or do_rowan):
        print("\n[dry run] nothing executed. Re-run with --cofold, --rowan or --all.")
        print("  --rowan spends Rowan credits; SQM dominates the cost.")
        return 0

    if do_cofold:
        print(f"\n=== co-folding ({a.poses}) -- GPU")
        rc = _run([sys.executable, str(ROOT / "scripts" / "run_predictors.py"), a.poses],
                  "needs CUDA_VISIBLE_DEVICES and the model env vars")
        if rc != 0:
            print("co-folding failed; not submitting poses that do not exist.")
            return rc

    if do_rowan:
        poses_exist = (WORK / "px_out").is_dir() or (WORK / "boltz_out").is_dir()
        if not poses_exist:
            print(f"\nno poses under {WORK.relative_to(ROOT)}/ -- the pose-based rescorers "
                  "(SQM/GNINA/AEV-PLIG) take a pose, not a sequence. Run with --cofold first "
                  "(or --all). _work/ is gitignored, so poses are regenerated, not fetched.\n"
                  "If you only want affinity columns, scripts/run_rowan_cofold.py folds on "
                  "Rowan and needs no GPU at all.")
            return 1
        cmd = [sys.executable, str(ROOT / "scripts" / "run_rowan.py"), "--submit",
               "--poses", a.poses]
        if a.max_credits is not None:
            cmd += ["--max-credits", str(a.max_credits)]
        else:
            print("\nrefusing to submit without --max-credits: an uncapped workflow can "
                  "spend without bound, and SQM is the expensive one.")
            return 1
        print("\n=== Rowan -- SPENDS CREDITS")
        rc = _run(cmd, "needs ROWAN_API_KEY")
        if rc != 0:
            return rc

    print("\n=== scoring")
    _run([sys.executable, str(ROOT / "scripts" / "score.py")])
    print("\nNext: commit predictions/ and push. The figures workflow re-scores and "
          "redraws on main, so the committed figure stays in step with the CSVs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
