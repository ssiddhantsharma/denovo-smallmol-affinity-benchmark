---
name: add-smallmol-paper
description: Add a published de-novo designed protein / small-molecule affinity paper to this benchmark. Runs extract -> verify -> deterministic append, then says what still needs a GPU.
---

# Add a paper to the de-novo small-molecule affinity benchmark

Pipeline: **extract (agent) → verify (agent) → append (deterministic Python) → predict (GPU, separate)**.

The agents never write `reference/`. Only `scripts/append_paper.py` does, which is what makes the set auditable: a hallucinated row has to survive checks that only look at values the extract already contains.

## Before starting

Is it in scope? A **designed protein** binding a **small molecule**, with a **measured** affinity. The commonest false positive is a paper designing a small-molecule *drug* against a natural protein — same words, opposite direction. If unsure, say so rather than spending a long extraction.

Already present? `cut -d, -f6 reference/system_reference.csv | sort -u` lists the source tags. Re-ingesting a paper wastes an extraction and risks a near-duplicate under a different id — `append_paper.py` dedups on (sequence, canonical SMILES), which catches the same pair but not a paraphrased construct name.

## Steps

1. **Extract.** Launch `smallmol-extractor` with the DOI/URL and a scratch output path. Relay its coverage summary — how many pairs it could not extract matters as much as how many it could.

2. **Verify.** Launch `smallmol-verifier` on the raw CSV, writing a cleaned one. Report kept / dropped / unverified. If everything dropped, stop and say so; do not ingest an empty file.

3. **Dry run.** `python scripts/append_paper.py --extract <cleaned>.csv --dry-run`. Read the rejections — they are the useful output. A "pKd outside range" is usually a unit slip, not a bad measurement.

4. **Append.** Same command without `--dry-run`. It backs both reference files up to `reference/backups/` first.

5. **Check.** `python scripts/qc_dataset.py`. It will now FAIL on "every predictions file covers exactly the reference ids" — that is correct and is the reminder for step 6, not a bug.

6. **Predict.** The new rows have no predictions, so they are in the reference set but in no figure. Run the predictors over them (`scripts/run_predictors.py` and friends). **This needs a GPU and is the only paid step** — confirm before launching it.

7. **Rescore and refresh.** `python scripts/score.py`, then `scripts/make_figure.py`. Re-run `qc_dataset.py`; it should pass now.

8. **Report.** New rows, new ligands, whether any ligand crossed 2 binders (which is what adds a within-ligand comparison), and the refreshed figure.

## Guardrails

- Never edit `reference/*.csv` by hand or let an agent write them.
- Re-running the same paper is safe: the dedup is on (sequence, canonical SMILES).
- A new ligand with only one binder adds an AUROC row but contributes nothing to the within-ligand correlation. Worth saying out loud in the report, because it looks like growth and partly is not.
