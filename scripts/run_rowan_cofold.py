"""Score the benchmark through Rowan's hosted co-folding: no GPU, no local poses.

`run_rowan.py` submits a pose this repo produced locally (Protenix/Boltz on a GPU),
split into apo protein + ligand. That path is still the right one for the
pose-based rescorers (SQM, GNINA, AEV-PLIG), which need a pose as input.

This path needs nothing local. Rowan's co-folding workflow takes a protein
SEQUENCE and a ligand SMILES, folds the complex server-side, and -- with
`ligand_binding_affinity_index` -- returns Boltz-2's affinity head alongside the
structure. So the whole chain for these columns is hosted, and the GPU box is
not in it.

What comes back per system (`CofoldingResult`):
    affinity_score.pred_value          Boltz-2 affinity, log10-ish; the repo's
                                       convention is pK = 6 - pred_value
    affinity_score.probability_binary  P(binder)
    scores.iptm / ptm / confidence     interface and global confidence
    pose, mmgbsa_score, strain, posebusters_valid

**These are NOT the repo's existing `boltz-2*` columns and must not be merged
into them.** Those were produced locally, at a pinned version, with this repo's
own MSA settings; Rowan's are a different deployment and a different MSA server.
Scoring one against the other as if they were one method is the provenance error
that makes a leaderboard quietly wrong, so every column written here is prefixed
`rowan-cofold-`.

It spends credits, so it is dry-run by default and every spending path is capped:

    python scripts/run_rowan_cofold.py                       # plan only
    python scripts/run_rowan_cofold.py --submit --max-credits 50 --max-systems 5
    python scripts/run_rowan_cofold.py --collect             # fetch finished workflows
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "reference" / "system_reference.csv"
PRED = ROOT / "predictions"
# TRACKED, not scratch. Two reasons. Submit and collect are separate runs because
# Rowan executes server-side, and a GitHub Actions runner is destroyed between
# them -- a state file under _work/ would lose every workflow id and the credits
# with it. And the map from a row to the workflow that produced its number is
# provenance worth keeping anyway.
STATE = ROOT / "reference" / "rowan_cofold_workflows.json"

# column -> how to pull it out of a CofoldingResult, and which way is better
COLUMNS = {
    "rowan-cofold-affinity": ("affinity", "higher = tighter (6 - pred_value)"),
    "rowan-cofold-pbind": ("pbind", "higher = more likely a binder"),
    "rowan-cofold-iptm": ("iptm", "higher = more confident interface"),
}


def rows() -> list[dict]:
    with REF.open() as fh:
        return list(csv.DictReader(fh))


def already_scored() -> set[str]:
    path = PRED / "rowan-cofold-affinity_predictions.csv"
    if not path.is_file():
        return set()
    with path.open() as fh:
        return {r["id"] for r in csv.DictReader(fh)}


def _state() -> dict:
    return json.loads(STATE.read_text()) if STATE.is_file() else {}


def _save_state(d: dict) -> None:
    STATE.parent.mkdir(parents=True, exist_ok=True)
    STATE.write_text(json.dumps(d, indent=2) + "\n")


def submit(todo: list[dict], model: str, max_credits: int, draft: bool) -> None:
    import rowan

    state = _state()
    for r in todo:
        if r["id"] in state:
            print(f"  already submitted: {r['id']} ({state[r['id']]})")
            continue
        wf = rowan.submit_protein_cofolding_workflow(
            initial_protein_sequences=[r["sequence"]],
            initial_smiles_list=[r["smiles"]],
            # index into initial_smiles_list: without it the workflow folds but
            # returns no affinity, which is the whole point of this path
            ligand_binding_affinity_index=0,
            model=model,
            use_msa_server=True,
            name=f"benchmark {r['id']}",
            max_credits=max_credits,
            is_draft=draft,
        )
        uuid = getattr(wf, "uuid", str(wf))
        print(f"  {'drafted' if draft else 'submitted'} {r['id']} -> {uuid}")
        if draft:
            # A draft executes nothing, so it has no result to collect. Recording
            # it would make the next REAL run skip this system as "already
            # submitted" and silently leave it unscored.
            continue
        state[r["id"]] = uuid
        _save_state(state)        # after each one: a crash must not lose the ids
    if draft:
        print("\n[draft] nothing recorded: a draft has no result to collect, and "
              "recording it would make the next real run skip these systems.")
    else:
        print(f"\n{len(state)} workflows tracked in {STATE.relative_to(ROOT)}")


def collect(wait: bool) -> None:
    import rowan

    state = _state()
    if not state:
        print(f"nothing submitted (no {STATE.relative_to(ROOT)})")
        return
    out: dict[str, dict[str, float]] = {c: {} for c in COLUMNS}
    pending, failed = [], {}
    for sid, uuid in state.items():
        try:
            wf = rowan.retrieve_workflow(uuid)
            if wait:
                wf.wait_for_result()
            if not getattr(wf, "is_finished", lambda: False)():
                pending.append(sid)
                continue
            res = wf.data if hasattr(wf, "data") else wf.result
            cof = res.cofolding_results[0] if getattr(res, "cofolding_results", None) else res
            aff = getattr(cof, "affinity_score", None)
            sco = getattr(cof, "scores", None)
            if aff and aff.pred_value is not None:
                out["rowan-cofold-affinity"][sid] = 6.0 - float(aff.pred_value)
            if aff and aff.probability_binary is not None:
                out["rowan-cofold-pbind"][sid] = float(aff.probability_binary)
            if sco and sco.iptm is not None:
                out["rowan-cofold-iptm"][sid] = float(sco.iptm)
        except Exception as exc:                      # noqa: BLE001 - one bad id must not lose the batch
            failed[sid] = f"{type(exc).__name__}: {exc}"

    for col, (_, note) in COLUMNS.items():
        vals = out[col]
        if not vals:
            print(f"  {col}: nothing returned yet")
            continue
        path = PRED / f"{col}_predictions.csv"
        with path.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["method", "id", "predicted_affinity"])
            for sid, v in sorted(vals.items()):
                w.writerow([col, sid, v])
        print(f"  wrote {path.relative_to(ROOT)}: {len(vals)} rows  ({note})")

    if pending:
        print(f"\nstill running ({len(pending)}): {', '.join(pending[:5])}"
              f"{' ...' if len(pending) > 5 else ''} -- re-run --collect later")
    for sid, why in failed.items():                  # surfaced, never silently dropped
        print(f"  FAILED {sid}: {why}")
    print("\nNote: scripts/qc_dataset.py requires every predictions file to cover every "
          "reference id. A partial file will fail that check until the rest return, "
          "which is the intended reminder.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--submit", action="store_true", help="SPENDS CREDITS")
    ap.add_argument("--collect", action="store_true", help="fetch finished workflows")
    ap.add_argument("--wait", action="store_true", help="with --collect, block until done")
    ap.add_argument("--model", default="boltz_2",
                    help="boltz_2 | boltz_2_1 | openfold_3 | chai_1r | esmfold_2")
    ap.add_argument("--max-credits", type=int, default=None, help="per-workflow cap")
    ap.add_argument("--max-systems", type=int, default=None, help="cap this session")
    ap.add_argument("--draft", action="store_true",
                    help="validate the full submission path without executing or spending")
    a = ap.parse_args()

    todo = [r for r in rows() if r["id"] not in already_scored()]
    print(f"{len(rows())} reference systems | {len(todo)} without a rowan-cofold score")
    if a.collect:
        collect(a.wait)
        return 0
    if not todo:
        print("nothing to submit.")
        return 0
    if a.max_systems is not None:
        todo = todo[: a.max_systems]
        print(f"capped to {len(todo)} this session")
    for r in todo[:8]:
        print(f"  {r['id']:42s} {len(r['sequence']):>4} aa  {r['smiles'][:40]}")
    if len(todo) > 8:
        print(f"  ... and {len(todo) - 8} more")

    if not a.submit:
        print("\n[dry run] nothing submitted. Add --submit with --max-credits.")
        print("  Each system is one hosted co-folding workflow; Boltz-2 with an MSA "
              "search is the cost driver.")
        return 0
    if a.max_credits is None and not a.draft:
        print("\nrefusing to submit without --max-credits. Rowan caps per workflow, so "
              "this and --max-systems are the only things bounding a session.")
        return 1
    if not os.environ.get("ROWAN_API_KEY"):
        print("\nROWAN_API_KEY is not set.")
        return 1
    print(f"\n=== submitting {len(todo)} co-folding workflows ({a.model})"
          f"{' [DRAFT: no execution, no credits]' if a.draft else ' -- SPENDS CREDITS'}")
    submit(todo, a.model, a.max_credits, a.draft)
    print("\nWorkflows run server-side; re-run with --collect when they finish.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
