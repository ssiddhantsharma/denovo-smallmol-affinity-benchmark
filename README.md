# de-novo small-molecule affinity benchmark

Do protein–ligand affinity predictors work on **de-novo designed** binders? Affinity leaderboards
score them on natural targets; this is a small complementary split built from published de-novo
designs.

**61 protein–ligand pairs:**
- **50 binders** with a measured KD (pKd).
- **11 matched negatives**: the same designed protein paired with a *wrong* ligand (confirmed
  non-binding), so a predictor has to read the interface, not the fold alone.

Each method is scored the same way for all 61 systems and joined to the ground truth by `id`
(cofolder affinity heads are converted to `pK = 6 - affinity_pred_value`).

## Results

![benchmark](figures/benchmark.png)

Left panel: rank affinity among the 50 binders (Spearman vs pKd). Right panel: tell the correct
ligand from the wrong one across all 61 (AUROC). Bars are 95% bootstrap CIs. This shows one metric
per model plus the baselines; the full leaderboard (all metrics) is in
[`figures/benchmark_full.png`](figures/benchmark_full.png), and `score.py` prints every metric.

- **Ranking affinity is hard.** The dedicated affinity heads do not beat lipophilicity: Boltz-2
  +0.30, Nesso-1 +0.27, versus cLogP +0.35. Structural confidence edges ahead (Boltz-2 ipTM +0.42);
  the interface-PAE metrics (min/mean ipae) cluster with it, none breaking out of the ~0.3–0.4 band.
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
non-binders from a handful of designed proteins. This small n reflects the field, not the search: a
systematic survey of the de novo literature turns up essentially no further designed binders of
*organic* small molecules with a precise measured KD — the remaining de novo binders are
metallo-cofactor systems (heme / Zn-porphyrin / Zn-chlorophyll maquettes) or designed pockets
grafted onto natural scaffolds. The set is therefore close to the available universe rather than a
sample of it. The three cofolders share the AF3-style family, so they are not independent. Cofolder
metrics carry run-to-run diffusion-sampling variation (~0.05 Spearman); values here are from a single
seeded fold. Positive KDs come from mixed assays and are not cross-calibrated.

## Methods scored

Every predictor is a third-party model or tool; this repo only wires them to the split, converts
their outputs to a common convention, and scores them. Please cite the original work when using any
result here.

| column(s) | method | family | source |
|---|---|---|---|
| `boltz-2`, `boltz-2-pbind`, `boltz-2-iptm`, `boltz-2-pae-{min,mean}` | Boltz-2 (affinity head + confidence) | co-folding | [jwohlwend/boltz](https://github.com/jwohlwend/boltz) |
| `nesso-1`, `nesso-1-pbind` | Nesso-1 (affinity head) | co-folding | [recursionpharma/nesso](https://github.com/recursionpharma/nesso) |
| `protenix-{iptm,ligiptm,gpde,ranking,pae-min,pae-mean}` | Protenix-v2 (confidence + ranking) | co-folding | [bytedance/Protenix](https://github.com/bytedance/Protenix) |
| `dtsfm-cosine` | dtSFM encoder cosine (drug–target specificity FM; MoLFormer-XL + ESM2 features) | sequence-native FM | see `scripts/run_dtsfm.py` |
| `clogp`, `molecular-weight` | Crippen cLogP, molecular weight | physicochemical baseline | [RDKit](https://www.rdkit.org) |
| *(planned, see below)* `rowan-sqm`, `rowan-gnina`, `rowan-aevplig` | SQM (PM6-D3H4X/COSMO2), GNINA CNN, AEV-PLIG | physics / docking / ML rescoring | [Rowan](https://docs.rowansci.com) |

Supporting libraries: structure parsing with [gemmi](https://gemmi.readthedocs.io) and
[RDKit](https://www.rdkit.org); the leave-one-out combination uses
[scikit-learn](https://scikit-learn.org) on [NumPy](https://numpy.org); the figures
(`benchmark.png`, `benchmark_full.png`) are drawn with [Matplotlib](https://matplotlib.org). The
Rowan rescoring below is driven through the [`rowan-python`](https://github.com/rowansci/rowan-python)
SDK; GNINA is [gnina/gnina](https://github.com/gnina/gnina) (McNutt et al., *J. Cheminform.* 2021),
and AEV-PLIG, NESSO and the SQM stack are used through Rowan (see their docs for primary references).
De-novo designs and labels are from the papers cited per row (`source`, `source_url`) in
`reference/system_reference.csv`.

## Rowan physics / docking rescoring (in progress)

The open question this split raises — *does any method actually track KD on de-novo binders, or does
that only work on natural complexes?* — is being probed with [Rowan](https://docs.rowansci.com)'s
binding-affinity workflow, which adds scorer families that are **not** co-folding heads:

- **SQM** — semi-empirical QM (PM6-D3H4X geometry optimization + COSMO2 single-point in water) on a
  truncated pocket. Physics-based, so not trained on natural protein–ligand data — the most direct
  test of whether the natural→de-novo shift is a *learning* artifact.
- **GNINA** — CNN docking affinity, and **AEV-PLIG** — an ML interaction-graph scorer. Both are
  trained on natural complexes (PDBbind-style), so a drop here would corroborate the same OOD wall
  from a different modeling family.

All three score a *bound pose*; the poses come from re-folding the 61 systems with Protenix-v2 (the
best-ranked of its 5 samples) and are fed to Rowan through the
[`rowan-python`](https://github.com/rowansci/rowan-python) SDK in its apo-protein + external-pose
mode. Wiring lives in `scripts/run_rowan.py`. Results and figures will be added once the runs
complete.

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

# Rowan physics/docking rescoring (needs Protenix poses in _work/px_out + a ROWAN_API_KEY)
python scripts/run_rowan.py                       # dry run: split poses, print the plan
ROWAN_API_KEY=... python scripts/run_rowan.py --submit --methods sqm,gnina,aevplig --max-credits 300
```

## Acknowledgements

The physics/docking rescoring is made possible by **compute credits generously provided by
[Rowan Scientific](https://rowansci.com)** — thank you. It also stands on the open models and tools
credited above (Boltz-2, Nesso-1, Protenix, dtSFM, RDKit, GNINA, AEV-PLIG, NESSO, and the
[`rowan-python`](https://github.com/rowansci/rowan-python) SDK); please cite the original work.
