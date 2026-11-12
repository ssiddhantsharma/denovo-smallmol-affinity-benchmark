# Contributing

Everything here is about *maintaining* the set. For what the benchmark finds, see the
[README](README.md).

## Adding a paper

The set is small because the field is: de-novo designed proteins with a *measured* affinity
for a small molecule are a few papers a year. The loop that grows it:

| | step | who | cost |
|---|---|---|---|
| 1 | `weekly` queries Europe PMC against a committed watermark, opens one ranked issue, and names the best new candidate | automatic, Monday | free |
| 2 | `ingest` extracts and verifies that candidate and **opens a PR** | automatic | API tokens |
| 3 | **review the PR** and merge | you | free |
| 4 | `weekly` collects last week's hosted runs and submits the new rows | automatic, Monday | ~46 credits |
| 4b | `gpu_autoscore.sh` on the GPU box scores the GPU-only columns | cron | GPU |
| 5 | `figures` re-scores, redraws, and restamps the README | automatic | free |

To add one by hand: run the `add-smallmol-paper` skill, or
`python scripts/ingest_paper.py --doi <doi> --out candidate.csv`, then
`python scripts/append_paper.py --extract candidate.csv --dry-run` before the real append.

**Step 3 is the one gate, and it stays.** `ingest` runs a model over the paper and writes a
candidate CSV. `append_paper.py` is still the only thing that writes `reference/`, and the job
**cannot push to main** — it pushes a branch and opens a PR that `qc` checks. A wrong row is
worse than a missing one, because nobody can tell which rows to distrust, so a person reads the
diff. The automation removes the fetching and transcribing, not the judgement.

**Expect most ingests to find nothing.** Europe PMC serves open-access full text; most of this
field's sequences sit in a PDF supplement behind a paywall. For a closed paper the extractor
correctly returns zero rows and opens no PR.

**The search query favours precision.** `"de novo" AND "small molecule"` is dominated by the
opposite problem — designing a *drug* against a natural protein — at about 2 useful hits in 23.
Requiring the designed thing to be the protein gives 5 hits over three years of which 2 are real;
adding biosensor terms for recall brought back SPR instruments and imprinted polymers. Recall is
the weak side on purpose: treat the weekly issue as a reading list.

## Scoring a new row

A merged row is in no figure until it is scored. `qc_dataset.py` reports the gap as a warning
naming the method and count; ids *not* in the reference fail instead, since that is always a bug.

**Hosted, no GPU.** Rowan's co-folding takes a sequence and a SMILES, folds server-side and
returns Boltz-2's affinity head with it:

```
gh secret set ROWAN_API_KEY                                     # once
gh workflow run rowan                                           # draft: validates, free
gh workflow run rowan -f draft=false -f max_systems=2 -f max_credits=50
gh workflow run rowan -f action=collect
```

Locally: `python scripts/run_rowan_cofold.py --submit --max-credits 50 --max-systems 5`.

Columns are `rowan-<model>-{affinity,pbind,iptm}`. **The model is in the name deliberately**:
`--model` is a dispatch input, so a fixed name would let a Chai-1r run overwrite Boltz-2's numbers
in place. They are also separate from the local `boltz-2*` columns and must stay separate — same
model family, different deployment and MSA server.

Measured cost is **9.1–9.4 credits per system**, so all 60 is roughly 560. Read it only from a
finished workflow: `credits_charged` grows while a run is in flight, and an early read of 3.06 is
how the first estimate came out three times too low.

**On the GPU box.** 11 of the 27 columns need one: the 9 `protenix-*` (Rowan hosts `boltz_2`,
`boltz_2_1`, `openfold_3`, `chai_1r`, `esmfold_2`, not Protenix) and `boltz-2-pae-{min,mean}` plus
ipDAE (these need the full residue-residue PAE matrix; a Rowan result carries `lddt`, `ptm`, `iptm`
and `confidence_score`, no PAE). Everything else runs hosted or on CPU.

Hosted runners have no GPU, and a *self-hosted* runner on a public repo would let any fork's pull
request execute code on that machine — so there is none. `scripts/gpu_autoscore.sh` inverts it: a
cron **on the GPU box** pulls main, runs `update.py --cofold`, runs QC and pushes. Nothing inbound
is exposed and no runner token exists.

```
# on the GPU box, crontab -e
30 10 * * 1  /path/to/repo/scripts/gpu_autoscore.sh >> ~/autoscore.log 2>&1
```

## Workflows

| workflow | trigger | does | cost |
|---|---|---|---|
| `ci` | push, PR | `ruff check scripts tests`, then `pytest` | free |
| `qc` | PR touching `reference/`, `predictions/`, `qc_dataset.py` | 21 integrity checks; row-count floor from the PR base branch | free |
| `figures` | push touching `predictions/`, `reference/`, `score.py`, `make_figure.py` | re-score, redraw, restamp, commit | free |
| `watch-sources` | `workflow_call`, dispatch | one Europe PMC query against a committed watermark; opens one issue. No model, no key, **cannot write `reference/`** | free |
| `ingest` | called by `weekly`, dispatch | extract + verify one paper and **open a PR**; cannot push to main; skips cleanly without `ANTHROPIC_API_KEY` | API tokens |
| `weekly` | Monday 09:00 UTC, dispatch | `watch-sources` → `ingest` the best candidate → collect finished Rowan runs → submit new rows | ~46 credits |
| `rowan` | dispatch only | submit or collect hosted co-folding, larger batches | credits |

Public repo, so Actions minutes are free; the only real cost is Rowan credits and ingest tokens.

**Collect runs before submit**, because a Rowan workflow finishes minutes to hours after
submission — the run that submits cannot also collect. Monday's collect turns last Monday's
submissions into predictions, and that push triggers `figures`. Without it the loop would submit
forever and never produce a number.

### What bounds the spending

Three things, because the first two do not bound a session:

- `--max-credits` caps one Rowan workflow. Rowan's own cap is per workflow, so this alone bounds
  nothing larger.
- `--max-systems` caps one run. Neither stops someone dispatching ten times.
- `BUDGET_CREDITS` in `run_rowan_cofold.py` is a ceiling on the **cumulative** charge. It reads the
  committed ledger in `reference/rowan_cofold_workflows.json` and refuses to submit past it. A
  recorded workflow whose charge has not been collected yet counts at the measured rate rather than
  as free, so an uncollected run is never invisible to it.

The Monday run spends up to 5 systems (~46 credits) plus ingest tokens; dispatch
`weekly -f score=false` to watch without scoring. At that rate the 800-credit ceiling is about four
months, after which submit refuses until someone raises it deliberately. Set **Max credit overage**
in the Rowan dashboard too: it defaults to unlimited, there is no billing API, and it is the only
thing that can bound spend this repo does not know about.

### Keys

Both live in repository secrets and are never printed — a workflow logs a key's length, not its
value. `rowan.yml` is `workflow_dispatch`-only, which GitHub restricts to write access. Fork pull
requests never receive repository secrets, and no workflow here uses `pull_request_target`, the
trigger that would hand a fork's code the keys. Dispatch inputs are passed through the environment
rather than interpolated into a shell command. `rowan.yml`'s `draft` defaults to **true**: the first
run of any change validates the whole submission path, executes nothing and costs nothing.

## Conventions

- `scripts/methods.py` is the registry: every column, what produces it, and its figure category.
  Adding a scorer means adding it there — a test and `qc_dataset.py` both fail on a registered
  method with no category, which is what stops `make_figure.py` dying on a bare `KeyError`.
- Prediction files are `predictions/<method>_predictions.csv` with `method,id,predicted_affinity`,
  joined to the reference by `id`. Lower-is-better metrics go in `LOWER_BETTER` in `score.py`.
- Never gate CI on a human-readable sentence. `run_rowan_cofold.py` and `update.py` both have
  `--count` for this: `weekly` used to grep a prose line, and renaming a column made the grep match
  nothing, read zero and skip the submit — a parse failure and "nothing to do" looked identical.
- `yaml.safe_load` is not workflow validation. A job block mis-indented into
  `on.workflow_dispatch.inputs` parses fine and GitHub rejects the file; so does a job requesting a
  permission its workflow lacks. Both surface as a run failure with no log, and both have tests.
- CI pins `matplotlib` to the version `figures.yml` uses, so the committed figure does not change
  for no reason anyone can see.
