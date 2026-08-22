# de-novo small-molecule affinity benchmark

Do protein-ligand affinity predictors work on **de-novo designed** binders? Existing leaderboards
(e.g. OpenBind EV-A71 2A) score them on natural targets. This is a small, clean, complementary split
built entirely from published de-novo designs: **47 protein-ligand pairs**, each a designed protein
plus a small molecule, with an experimental label.

- **37 cognate binders** with a measured KD (pKd), from published de-novo campaigns.
- **10 matched-specificity negatives**: the same designed proteins paired with a *wrong* ligand
  (confirmed non-binding). This makes the negatives hard: the protein is identical, only the ligand
  is wrong, so a predictor must read the interface, not the fold.

Six predictor families are scored the same way for both cofolders, `pK = 6 - affinity_pred_value`
(the OpenBind convention): **Boltz-2** (affinity head + iptm), **Nesso-1**, **Protenix-v2** (iptm,
gpde, ranking), and the trivial **MW / cLogP** baselines any real signal must beat.

## Two tasks

![benchmark](figures/benchmark.png)

*Rank affinity magnitude among the 37 binders (left, Spearman rho vs pKd) and tell the cognate ligand
from the wrong one on all 47 (right, AUROC). Bars are 95% bootstrap CIs; dashed lines are the cLogP
baseline (left) and chance (right); the star is a leave-one-out combination of four metrics.*

**Ranking affinity magnitude is hard.** The dedicated affinity heads do not beat lipophilicity:
Boltz-2 pK rho = +0.39, Nesso-1 pK = +0.33, versus cLogP = +0.40. Structural proxies edge it
(Protenix / Boltz iptm ~ +0.48). Every CI is wide (n=37), so no single metric is clearly separated
from cLogP. RMSE on pK is 1.3-1.5 here, versus Nesso's 0.86 on the natural EV-A71 target: a
distribution-shift gap.

**Discriminating binding is easy.** The cofolders cleanly separate the cognate ligand from the wrong
one (Boltz-2 iptm AUROC = 0.94, Protenix iptm / ranking = 0.89, Boltz pbind = 0.87), while the
baselines are at or below chance (cLogP 0.44, MW 0.24). So the models genuinely read the interface,
they just cannot rank affinity *magnitude* among true binders.

**Metrics combine.** A leave-one-out combination of Boltz-2 pK, Nesso-1 pK, Protenix lig-iptm, and
cLogP reaches rho = +0.59 for ranking, above any single metric (specificity does not improve over
Boltz-2 iptm alone). The signals are complementary.

| task | best single | cLogP baseline | combination (LOO) |
|---|---|---|---|
| rank pKd (37, Spearman) | +0.48 (protenix iptm) | +0.40 | **+0.59** |
| specificity (47, AUROC) | 0.94 (boltz iptm) | 0.44 | 0.83 |

## Caveats

- **Small n** (37 binders, 10 negatives). CIs are wide; treat single-metric orderings as suggestive.
- **Related judges.** Boltz-2, Nesso-1, and Protenix share the AF3-style cofolding family; they are
  not independent oracles.
- **Curated panel**, not an unbiased sample: it is the de-novo designs with public sequences,
  ligands, and labels that could be assembled.
- **In silico labels for negatives** are experimental non-binders from the source papers, but the
  positives' KDs come from several assays (ITC / SPR / MST / fluorescence) and are not cross-calibrated.

## Reproduce

```
# 1. run each predictor (cofolders need their own GPU env) + the RDKit baselines
CUDA_VISIBLE_DEVICES=0 BOLTZ=/path/to/boltz            python scripts/run_predictors.py boltz
CUDA_VISIBLE_DEVICES=0 NESSO=/path/to/nesso            python scripts/run_predictors.py nesso
CUDA_VISIBLE_DEVICES=0 PROTENIX_DIR=/path/to/Protenix  python scripts/run_predictors.py protenix
python scripts/run_predictors.py baselines
# 2. score + figure
python scripts/score.py
python scripts/make_figure.py
```

Layout follows the OpenBind EV-A71 2A benchmark:

```
reference/experimental_reference_ground_truth.csv   id, is_binder, experimental_pKD
reference/system_reference.csv                       id, protein, ligand_name, sequence, smiles, source, source_url
predictions/<method>_predictions.csv                 method, id, predicted_affinity   (one file per method, strict 3 columns)
```

Predictions join to the reference by `id`; every `predictions/*.csv` is scored automatically.

## Credit

Predictors: [Boltz-2](https://github.com/jwohlwend/boltz), [Nesso-1](https://github.com/recursionpharma/nesso),
[Protenix](https://github.com/bytedance/Protenix). Benchmark design follows the
[OpenBind Consortium](https://github.com/OpenBind-Consortium) EV-A71 2A leaderboard. De-novo designs
and labels are from the source papers cited per row in `data/bench47.json`.
