# de-novo small-molecule affinity benchmark

[![ci](https://github.com/ssiddhantsharma/denovo-smallmol-affinity-benchmark/actions/workflows/ci.yml/badge.svg)](https://github.com/ssiddhantsharma/denovo-smallmol-affinity-benchmark/actions/workflows/ci.yml)

Do protein-ligand affinity predictors work on **de-novo designed** binders? Affinity leaderboards
score them on natural targets; this is a small complementary split built from published de-novo
designs.

**61 protein-ligand pairs:**
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
  the interface-PAE metrics (min/mean ipae) cluster with it, none breaking out of the 0.3 to 0.4
  band. CIs are wide at n=50.
- **Combining does not help.** A leave-one-out combination of four metrics (+0.37) is statistically
  tied with the best single metric, Boltz-2 ipTM (+0.42): the bootstrap CI on the difference is
  [-0.31, +0.23], spanning 0. For specificity it trails Boltz-2 ipTM (0.84 vs 0.91). No evidence the
  metrics carry complementary signal here.
- **A sequence-native affinity FM fails too; the failure is the data regime, not co-folding.**
  dtSFM (`dtsfm-cosine`), a 714k-pair drug-target specificity foundation model, drops to chance on de
  novo: AUROC 0.53 (right vs wrong ligand) and Spearman +0.19 (pKd), *below* the cLogP baseline. Yet
  its featurization validates on in-distribution natural pairs (true-vs-shuffled AUROC 0.863; see
  `scripts/run_dtsfm.py`, which gates on that positive control before scoring). So neither a
  structural co-folder nor a sequence-native FM escapes the natural-to-de-novo distribution shift.

| task | best single | cLogP | combination (LOO) |
|---|---|---|---|
| rank pKd (Spearman) | +0.42 | +0.35 | +0.37 (tied) |
| tell right from wrong ligand (AUROC) | 0.91 | 0.54 | 0.84 |

## Caveats

Small n (50 + 11): CIs are wide and orderings are suggestive. The 11 negatives are confirmed
non-binders from a handful of designed proteins. This small n reflects the field, not the search: a
systematic survey of the de novo literature turns up essentially no further designed binders of
*organic* small molecules with a precise measured KD. The remaining de novo binders are
metallo-cofactor systems (heme, Zn-porphyrin, Zn-chlorophyll maquettes) or designed pockets grafted
onto natural scaffolds, so the set is close to the available universe rather than a sample of it. The
three cofolders share the AF3-style family, so they are not independent. Cofolder metrics carry
run-to-run diffusion-sampling variation (~0.05 Spearman); values here are from a single seeded fold.
Positive KDs come from mixed assays and are not cross-calibrated.

## Methods scored

Every predictor is a third-party model or tool; this repo only wires it to the split, converts its
outputs to a common convention, and scores it. Full citations are in [References](#references);
please cite the original work when using any result here.

| column(s) | method | family | reference |
|---|---|---|---|
| `boltz-2`, `boltz-2-pbind`, `boltz-2-iptm`, `boltz-2-pae-{min,mean}` | Boltz-2 (affinity head + confidence) | co-folding | Passaro et al. 2025 |
| `nesso-1`, `nesso-1-pbind` | Nesso-1 (affinity head) | co-folding | Recursion / Valence Labs |
| `protenix-{iptm,ligiptm,gpde,ranking,pae-min,pae-mean}` | Protenix v2 (confidence + ranking) | co-folding | ByteDance AML 2025 |
| `dtsfm-cosine` | dtSFM encoder cosine (drug-target specificity FM) | sequence-native FM | dtSFM-v3, BIIE ETH Zürich |
| `clogp`, `molecular-weight` | Crippen cLogP, molecular weight | physicochemical baseline | RDKit |
| `rowan-sqm`, `rowan-gnina`, `rowan-aevplig` *(in progress)* | SQM (PM6-D3H4X/COSMO2), GNINA, AEV-PLIG | physics / docking / ML rescoring | via Rowan |

De-novo designs and labels are from the papers cited per row (`source`, `source_url`) in
`reference/system_reference.csv`.

## Rowan physics / docking rescoring (in progress)

The open question this split raises, *does any method actually track KD on de-novo binders, or does
that only work on natural complexes?*, is being probed with [Rowan](https://docs.rowansci.com)'s
binding-affinity workflow, which adds scorer families that are not co-folding heads:

- **SQM**: semi-empirical QM (PM6-D3H4X geometry optimization, COSMO2 single-point in water) on a
  truncated pocket. Physics-based, so not trained on natural protein-ligand data, and the most direct
  test of whether the natural-to-de-novo shift is a *learning* artifact.
- **GNINA** (CNN docking affinity) and **AEV-PLIG** (an ML interaction-graph scorer). Both are
  trained on natural complexes (PDBbind-style), so a drop here would corroborate the same
  out-of-distribution wall from a different modeling family.

All three score a *bound pose*. The poses come from re-folding the 61 systems with Protenix v2 (the
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

## Tests

```
pip install numpy scikit-learn pytest
pytest -q
```

`tests/` covers the scoring statistics (Spearman, Pearson, AUROC, RMSE), dataset integrity
(61 systems = 50 binders + 11 negatives, every prediction joins to the ground truth), and that
`score.py` runs end to end. CI (ruff lint + pytest) runs on every push.

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

## References

Models and scorers:

- **Boltz-2**: Passaro, Corso, Wohlwend, et al. *Boltz-2: Towards Accurate and Efficient Binding
  Affinity Prediction.* bioRxiv (2025). Code: [jwohlwend/boltz](https://github.com/jwohlwend/boltz).
- **Nesso-1**: Recursion / Valence Labs. *Nesso-1* technical report and model card:
  [recursionpharma/nesso](https://huggingface.co/recursionpharma/nesso).
- **Protenix v2**: ByteDance AML AI4Science Team (Chen et al.). *Protenix: Advancing Structure
  Prediction Through a Comprehensive AlphaFold3 Reproduction.* bioRxiv (2025),
  [doi:10.1101/2025.01.08.631967](https://doi.org/10.1101/2025.01.08.631967).
- **dtSFM**: Specificity Foundation Model, BIIE, ETH Zürich. Weights
  [SFM-BIIE-ETHZ/dtSFM-v3](https://huggingface.co/SFM-BIIE-ETHZ/dtSFM-v3); features from MoLFormer-XL
  (Ross et al., *Nat. Mach. Intell.* 2022) and ESM2 (Lin et al., *Science* 2023).
- **GNINA**: McNutt, Francoeur, Aggarwal, et al. *GNINA 1.0: molecular docking with deep learning.*
  *J. Cheminform.* 13, 43 (2021). Code: [gnina/gnina](https://github.com/gnina/gnina).
- **AEV-PLIG**: Valsson, Warren, Deane, Magarkar, Morris, Biggin. *Narrowing the gap between machine
  learning scoring functions and free energy perturbation using augmented data.* (2025),
  [PMC11807228](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11807228/).
- **SQM (PM6)**: Stewart. *Optimization of parameters for semiempirical methods V.* *J. Mol. Model.*
  13, 1173 (2007). Run via Rowan as PM6-D3H4X with COSMO2 solvation.

Software:

- **rowan-python** SDK: Rowan Scientific,
  [rowansci/rowan-python](https://github.com/rowansci/rowan-python).
- **RDKit**: Landrum et al. *RDKit: Open-source cheminformatics*, [rdkit.org](https://www.rdkit.org).
- **gemmi**: Wojdyr. *GEMMI: A library for structural biology.* *J. Open Source Softw.* 7, 4200 (2022).
- **scikit-learn**: Pedregosa et al. *Scikit-learn: Machine Learning in Python.* *JMLR* 12, 2825 (2011).
- **NumPy**: Harris et al. *Array programming with NumPy.* *Nature* 585, 357 (2020).
- **Matplotlib**: Hunter. *Matplotlib: A 2D graphics environment.* *Comput. Sci. Eng.* 9, 90 (2007).

## Acknowledgements

The physics/docking rescoring is made possible by compute credits generously provided by
[Rowan Scientific](https://rowansci.com); thank you. See [References](#references) for the models,
scorers, and libraries this benchmark builds on.
