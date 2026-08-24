# de-novo small-molecule affinity benchmark

Do protein-ligand affinity predictors work on **de-novo designed** binders? Affinity leaderboards
score them on natural targets; this is a small complementary split built from published de-novo
designs.

**47 protein-ligand pairs:**
- **37 binders** with a measured KD (pKd).
- **10 matched negatives**: the same designed protein paired with a *wrong* ligand (confirmed
  non-binding), so a predictor has to read the interface, not just the fold.

Six predictors, scored two ways (cofolders converted to `pK = 6 - affinity_pred_value`):
[Boltz-2](https://github.com/jwohlwend/boltz) (affinity head + iptm),
[Nesso-1](https://github.com/recursionpharma/nesso),
[Protenix-v2](https://github.com/bytedance/Protenix) (iptm, gpde, ranking), and MW / cLogP baselines.

## Results

![benchmark](figures/benchmark.png)

Left panel: rank affinity among the 37 binders (Spearman vs pKd). Right panel: tell the correct
ligand from the wrong one across all 47 (AUROC). Bars are 95% bootstrap CIs. The figure shows one
metric per model plus the baselines; `score.py` reports all eleven.

- **Ranking affinity is hard.** The dedicated affinity heads do not beat lipophilicity: Boltz-2
  +0.39, Nesso-1 +0.33, versus cLogP +0.40. Structural confidence edges ahead (iptm ~+0.48). CIs are
  wide at n=37.
- **Discriminating binding is easy.** Cofolders separate the correct ligand from the wrong one
  cleanly (Boltz-2 iptm AUROC 0.94), while the baselines sit at chance.
- **Metrics combine.** A leave-one-out combination reaches +0.59 for ranking, above any single
  metric; specificity does not improve over Boltz-2 iptm alone.

| task | best single | cLogP | combination (LOO) |
|---|---|---|---|
| rank pKd (Spearman) | +0.48 | +0.40 | **+0.59** |
| tell right from wrong ligand (AUROC) | 0.94 | 0.44 | 0.83 |

## Caveats

Small n (37 + 10): CIs are wide and orderings are suggestive. The 10 negatives are confirmed
non-binders from ~6 designed proteins. The combination's edge over the best single metric is modest.
The three cofolders share the AF3-style family, so they are not independent. Positive KDs come from
mixed assays and are not cross-calibrated.

## Layout

```
reference/experimental_reference_ground_truth.csv   id, is_binder, experimental_pKD
reference/system_reference.csv                       id, protein, ligand_name, sequence, smiles, source, source_url
predictions/<method>_predictions.csv                 method, id, predicted_affinity
```

Every `predictions/*.csv` joins to the reference by `id` and is scored automatically.

## Reproduce

```
# each cofolder needs its own GPU env; baselines are RDKit
CUDA_VISIBLE_DEVICES=0 BOLTZ=/path/to/boltz            python scripts/run_predictors.py boltz
CUDA_VISIBLE_DEVICES=0 NESSO=/path/to/nesso            python scripts/run_predictors.py nesso
CUDA_VISIBLE_DEVICES=0 PROTENIX_DIR=/path/to/Protenix  python scripts/run_predictors.py protenix
python scripts/run_predictors.py baselines
python scripts/score.py
python scripts/make_figure.py
```

## Credit

De-novo designs and labels are from the papers cited per row (`source`, `source_url`) in
`reference/system_reference.csv`.
