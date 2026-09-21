# de-novo small-molecule affinity benchmark

Do protein-ligand affinity predictors work on **de-novo designed** binders? Affinity leaderboards
score them on natural targets; this is a small complementary split built from published de-novo
designs.

**61 protein-ligand pairs:**
- **50 binders** with a measured KD (pKd).
- **11 matched negatives**: the same designed protein paired with a *wrong* ligand (confirmed
  non-binding), so a predictor has to read the interface, not the fold alone.

Five model families, scored two ways (cofolder affinity heads converted to `pK = 6 -
affinity_pred_value`): [Boltz-2](https://github.com/jwohlwend/boltz) (affinity head, ipTM, interface
PAE), [Nesso-1](https://github.com/recursionpharma/nesso) (affinity head),
[Protenix-v2](https://github.com/bytedance/Protenix) (ipTM, gPDE, ranking, interface PAE), and MW /
cLogP baselines, for 15 metric columns in all.

## Results

![benchmark](figures/benchmark.png)

Left panel: rank affinity among the 50 binders (Spearman vs pKd). Right panel: tell the correct
ligand from the wrong one across all 61 (AUROC). Bars are 95% bootstrap CIs. This shows one metric
per model plus the baselines; the full leaderboard (all 15 metrics) is in
[`figures/benchmark_full.png`](figures/benchmark_full.png), and `score.py` prints every metric.

- **Ranking affinity is hard.** The dedicated affinity heads do not beat lipophilicity: Boltz-2
  +0.30, Nesso-1 +0.27, versus cLogP +0.35. Structural confidence edges ahead (Boltz-2 ipTM +0.42);
  the interface-PAE metrics (min/mean ipae) cluster with it, none breaking out of the ~0.3-0.4 band.
  CIs are wide at n=50.
- **Combining does not help.** A leave-one-out combination of four metrics (+0.37) is statistically
  tied with the best single metric, Boltz-2 ipTM (+0.42): the bootstrap CI on the difference is
  [-0.31, +0.23], spanning 0. For specificity it trails Boltz-2 ipTM (0.84 vs 0.91). No evidence the
  metrics carry complementary signal here.
- **A sequence-native affinity FM fails too — the failure is the data regime, not co-folding.**
  dtSFM (`dtsfm-cosine`), a 714k-pair drug–target specificity foundation model, drops to chance on de
  novo: AUROC 0.53 (right vs wrong ligand) and Spearman +0.19 (pKd), *below* the cLogP baseline —
  even though its featurization validates on in-distribution natural pairs (true-vs-shuffled AUROC
  0.863; see `scripts/run_dtsfm.py`, which gates on that positive control before scoring). So neither a
  structural co-folder nor a sequence-native FM escapes the natural→de-novo distribution shift.

| task | best single | cLogP | combination (LOO) |
|---|---|---|---|
| rank pKd (Spearman) | +0.42 | +0.35 | +0.37 (tied) |
| tell right from wrong ligand (AUROC) | 0.91 | 0.54 | 0.84 |

## Caveats

Small n (50 + 11): CIs are wide and orderings are suggestive. The 11 negatives are confirmed
non-binders from a handful of designed proteins. This small n reflects the field, not the search: a systematic
survey of the de novo literature turns up essentially no further designed binders of *organic* small
molecules with a precise measured KD — the remaining de novo binders are metallo-cofactor systems (heme /
Zn-porphyrin / Zn-chlorophyll maquettes) or designed pockets grafted onto natural scaffolds. The set is
therefore close to the available universe rather than a sample of it. The combination does not beat the best single metric (tied
within bootstrap noise). The three cofolders share the AF3-style family, so they are not independent. Cofolder metrics carry
run-to-run diffusion-sampling variation (~0.05 Spearman); values here are from a single seeded fold.
Positive KDs come from mixed assays and are not cross-calibrated.

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
