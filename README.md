# de-novo small-molecule affinity benchmark

[![ci](https://github.com/ssiddhantsharma/denovo-smallmol-affinity-benchmark/actions/workflows/ci.yml/badge.svg)](https://github.com/ssiddhantsharma/denovo-smallmol-affinity-benchmark/actions/workflows/ci.yml)

Do protein-ligand affinity predictors work on **de-novo designed** binders? Affinity leaderboards
score them on natural targets; this is a small complementary split built from published de-novo
designs.

**60 protein-ligand pairs:**
- **49 binders** with a measured KD (pKd).
- **11 confirmed non-binders**: 8 are a binder's protein paired with a *wrong* ligand, 3 are
  non-binding design variants. Either way a predictor has to read the interface, not the fold alone.

Affinity is ranked **within-ligand**: among the de-novo proteins that bind the *same* ligand, Spearman
of a metric against pKd, median across ligands (the analog of the within-target correlation used on
de-novo protein sets). It is reported in place of a pooled Spearman across all binders because the
set spans 19 ligands, and a pooled correlation is dominated by ligand composition, which cLogP
measures, rather than by the interface. A ligand-only score such as cLogP is
constant within a ligand and drops out of this statistic. Classification (right vs wrong ligand) is
one AUROC over all 60. Cofolder affinity heads are converted to `pK = 6 - affinity_pred_value`.

## Results

![benchmark](figures/benchmark.png)

<!-- figure-stamp -->
*Figure regenerated **2026-11-12 12:16 UTC** from 60 systems, 25 methods drawn.*
<br>*Not drawn (no AUROC until both classes are scored): `rowan-boltz-2-affinity`, `rowan-boltz-2-iptm`, `rowan-boltz-2-pbind`.*
<!-- /figure-stamp -->

Every metric, two panels: within-ligand rank-pKd (left, composition controlled; the dashed line is no
correlation, rho = 0) and right-vs-wrong-ligand AUROC (right; the dashed line is chance, 0.5). Bars
are 95% bootstrap CIs. `score.py` prints every value, including the pooled Spearman, which is
retained for transparency.

- **Co-folder interface confidence ranks highest within ligand.** Holding the molecule fixed
  (median over the 5 ligands with at least 3 protein binders), the interface-confidence metrics
  lead: Protenix ipTM and lig-ipTM at +0.70, Boltz-2 ipTM at +0.60, Protenix pLDDT-ligand at +0.43,
  with interface-PAE in the same band. The physics and machine-learning rescorers are mid-ranked
  (SQM +0.33, AEV-PLIG +0.35, GNINA +0.30), and the affinity heads rank below them (Boltz-2 +0.30,
  Nesso-1 -0.05). This is the same confidence family (iPTM, pTM, pLDDT, interface-PAE) that ranks
  KD on de-novo *protein* sets. Confidence intervals are wide, over 5 ligand groups of 3 to 9
  binders each, so the ordering is suggestive rather than decisive.
- **Almost none of that ordering is statistically separable.** `score.py` runs a paired bootstrap
  over the ligands that *both* metrics scored, and **18 of 21** metrics are not distinguishable
  from the leading one. The three it separates from with an interval clear of zero are Nesso-1
  P(bind), SQM and Nesso-1, and that is before correcting for 21 comparisons.
  The left panel supports an ordering of method families rather than of individual metrics: the
  gap between Protenix ipTM at +0.70 and Boltz-2 ipTM at +0.60 is two medians over the same five
  ligands.
- **The co-folders also lead on specificity.** On distinguishing the correct ligand from the wrong
  one, the interface-confidence metrics rank highest: Boltz-2 and Protenix ipTM at 0.85 to 0.91,
  Protenix pLDDT-ligand at 0.86, pTM at 0.79. cLogP, dtSFM and the physics and machine-learning
  rescorers are near chance. A pure-geometry pose baseline, ligand burial
  (`interface-ligburial`, the fraction of the ligand within 4.5 A of the protein), reaches 0.71,
  above cLogP, but carries no affinity-ranking signal that survives its interval (within-ligand rho
  +0.32 over 3 ligands, CI [-0.62, +0.71]). The correct ligand is more buried, but burial does not
  demonstrably order KD.
- **A sequence-native affinity model also fails, which implicates the data regime rather than
  co-folding.** dtSFM (`dtsfm-cosine`), a 714k-pair drug-target specificity foundation model, is at
  chance on de novo designs (AUROC 0.53, within-ligand +0.10), while its featurization validates on
  in-distribution natural pairs (true-versus-shuffled AUROC 0.863). `scripts/run_dtsfm.py` gates on
  that positive control first.

| task | best single | cLogP | combination (LOO) |
|---|---|---|---|
| rank pKd, within-ligand (median Spearman) | +0.70 | n/a (constant per ligand) | not defined |
| tell right from wrong ligand (AUROC) | 0.91 | 0.53 | 0.84 |

> **Correction, 2026-11-11.** The rank statistics broke ties by row order, using a Spearman
> shortcut that is not valid in the presence of ties. Five metrics carry ties, so their values
> were incorrect; the largest change is `interface-ligburial`'s within-ligand rho, published as
> +0.10 and corrected to +0.32. Ranks are now midranks, agreeing with scipy to 1e-16. Nine of 68
> point estimates changed and no conclusion was affected.

## Limitations

**The within-ligand axis rests on less data than the 60 pairs suggest.** Only 5 of the 19 ligands
have the 3 or more protein binders a within-ligand rank requires, so that median is taken over
5 ligands and 30 pairs for most metrics, and over 3 ligands and 18 pairs for `interface-ligburial`,
for which some poses are missing. The figure prints this n beside every bar and `score.py` reports
it as `n_lig` and `n_pairs`. The remaining 14 ligands contribute to the AUROC panel, which uses all
60 pairs, but not to the affinity panel.

**The set is small and the intervals are correspondingly wide.** The 11 negatives are confirmed
non-binders drawn from a small number of designed proteins. The size reflects the field rather than
the search: a systematic survey of the de novo literature finds few further designed binders of
*organic* small molecules with a precise measured KD. The remainder are metallo-cofactor systems
(heme, Zn-porphyrin and Zn-chlorophyll maquettes) or designed pockets grafted onto natural
scaffolds, so the set approximates the available universe rather than sampling it.

**The predictors are not independent, and the values are not stable.** The three cofolders share
the AF3-style architecture. Cofolder metrics carry run-to-run diffusion-sampling variation, which
on this set was measured rather than estimated and reaches 0.26 in within-ligand Spearman (see
[Sensitivity to re-folding](#sensitivity-to-re-folding)). Positive KDs come from mixed assays and
are not cross-calibrated.

## Methods scored

Every predictor is a third-party model or tool. This repository connects it to the split, converts
its outputs to a common convention and scores it. Full citations are in
[References](#references); please cite the original work when using any result here.

| column(s) | method | family | reference |
|---|---|---|---|
| `boltz-2`, `boltz-2-pbind`, `boltz-2-iptm`, `boltz-2-pae-{min,mean}` | Boltz-2 (affinity head + confidence) | co-folding | Passaro et al. 2025 |
| `nesso-1`, `nesso-1-pbind` | Nesso-1 (affinity head) | co-folding | Recursion / Valence Labs |
| `protenix-{iptm,ligiptm,ptm,plddt-binder,plddt-ligand,gpde,ranking,pae-min,pae-mean}` | Protenix v2 (confidence + ranking) | co-folding | ByteDance AML 2025 |
| `dtsfm-cosine` | dtSFM encoder cosine (drug-target specificity FM) | sequence-native FM | dtSFM-v3, BIIE ETH Zürich |
| `clogp`, `molecular-weight` | Crippen cLogP, molecular weight | physicochemical baseline | RDKit |
| `rowan-sqm`, `rowan-gnina`, `rowan-aevplig` | SQM (PM6-D3H4X/COSMO2), GNINA, AEV-PLIG | physics / docking / ML rescoring | via Rowan |
| `rowan-boltz-2-{affinity,pbind,iptm}` | Boltz-2 hosted end to end (sequence + SMILES in, no GPU) | co-folding | via Rowan |
| `interface-ncontacts`, `interface-ligburial` | protein-ligand contacts and ligand burial from the pose | interface geometry | `scripts/run_interface.py` (gemmi) |

The confidence metrics cover the de-novo metric set used on protein binders (iPTM, pTM, pLDDT,
interface-PAE, contacts). ipSAE, ipDAE and LIS, the strongest there, are protein-protein scores
built from the full residue-residue PAE matrix. Availability depends on the folder:

| folder | full residue-residue PAE | consequence |
|---|---|---|
| Boltz-2 | saved (`--write_full_pae`) | ipDAE is computable without a GPU, via `scripts/run_ipdae.py` |
| Protenix | not saved, only `chain_pair_pae_{min,mean}` | ipSAE, ipDAE and LIS would require a re-fold |

For a few-atom ligand the natural analog of these scores is lig-ipTM plus interface-PAE, which are
included for both folders.

De-novo designs and labels are from the papers cited per row (`source`, `source_url`) in
`reference/system_reference.csv`.

## Sensitivity to re-folding

The Protenix predictions were re-derived from the folds that still exist (`_work/px_out`, seed 101,
5 diffusion samples per system, averaged). They did not reproduce what was committed. The committed
values came from an earlier fold that no longer exists: regenerating from the surviving outputs
matches the current files to 0.000000 and the earlier release only to 0.036.

The earlier set also mixed two folds. `protenix-ptm`, both pLDDT columns and both
interface-geometry columns regenerated byte-identical, having been derived later from the surviving
fold, while `protenix-{iptm,ligiptm,gpde,ranking,pae-*}` did not. The figure therefore compared
metrics from one fold against metrics from another.

All Protenix-derived columns are now regenerated from the single surviving fold via
`python scripts/run_predictors.py protenix-parse`, which requires no GPU. The effect on the
reported values:

| metric | earlier release | surviving fold | delta |
|---|---|---|---|
| `protenix-gpde` | +0.50 | **+0.24** | **-0.26** |
| `protenix-ranking` | +0.47 | +0.58 | +0.11 |
| `protenix-pae-min` | +0.52 | +0.45 | -0.07 |
| `protenix-iptm` / `-ligiptm` | +0.65 | +0.70 | +0.05 |

Re-running the same model at the same seed moves a within-ligand median by up to 0.26, which
exceeds almost every gap in the ordering such values are used to justify. This is independent
support for the separability result above: where one re-fold moves a metric by 0.26, an ordering
whose neighbours differ by 0.03 to 0.15 is not resolved. Single-fold affinity benchmark values,
here and elsewhere, should be read as carrying uncertainty of this magnitude whether or not it is
drawn. The leading metric is unchanged (`protenix-iptm`), 18 of 21 metrics remain non-separable
from it, and the AUROC panel is largely unaffected.

The Boltz-2 and Nesso-1 columns are untouched: their outputs are gone, so they could not be
re-derived and were not re-run. Those columns are still from their original single fold.

## ipDAE

`scripts/run_ipdae.py` ports BindCraft2's `i_pDAE` to a protein-ligand interface. The
implementation is unit-tested but has not been run against the folds, which remain on the machine
that produced them, so no result above includes it.

Its contribution is expected to be small. The partner-count scale `d0` exceeds 1.0 only at 27 or
more partners, verified numerically and pinned by a test, and a ligand atom does not contact 27
residues. `d0` therefore pins to 1.0 and the score reduces to the maximum over anchors of mean
`1/(1+PAE²)`. This is not a monotone transform of mean interface PAE, but it reads the same matrix
as `boltz-2-pae-{min,mean}`, already reported at within-ligand +0.62 and +0.50. Given that 18 of 21
metrics are not separable from the leading one, ipDAE is more likely to fall within that group than
to add an independent axis. `--report` prints the measured fraction of anchors with `d0 > 1`.

Contact is defined as heavy-atom within 4.5 Å rather than CA-CA within 8 Å, since a ligand has no
CA, matching the definition `run_interface.py` uses. The script exits with an error if the
structure's residue count disagrees with the sequence, or if the PAE token count is not
`len(protein) + n_ligand_heavy_atoms`. Boltz-2's PAE is paired with Boltz-2's own structure, and
combining it with the Protenix pose would score one model's confidence against another model's
coordinates.

## Physics and docking rescoring

Whether any method tracks KD on de-novo binders, or only on natural complexes, is examined here
with [Rowan](https://docs.rowansci.com)'s binding-affinity workflow across scorer families that are
not co-folding heads. Each scores a bound pose, the best-ranked of Protenix v2's five samples per
system, submitted through the [`rowan-python`](https://github.com/rowansci/rowan-python) SDK in its
apo-protein plus external-pose mode. SQM additionally runs a protein-preparation step first.
Implementation is in `scripts/run_rowan.py`.

- **Semiempirical physics is mid-ranked once the ligand is controlled.** Within-ligand, SQM reaches
  +0.33, below the interface-confidence metrics at approximately +0.70. On the confounded pooled
  Spearman it reaches +0.40. That margin comes from ligand composition, the same property cLogP
  measures, rather than from the interface.
- **The learned rescorers are also mid-ranked.** GNINA at +0.30 and AEV-PLIG at +0.35 within-ligand,
  both trained on natural complexes of the PDBbind type.
- **None of the three contributes to specificity.** On distinguishing the correct ligand from the
  wrong one, SQM reaches 0.57, GNINA 0.58 and AEV-PLIG 0.47, all near chance and well below
  co-folder ipTM at 0.91.

On this set, physics and docking rescoring do not improve on the co-folders' own
interface-confidence metrics.

Limitations: SQM's truncated-pocket energies are noisy (cut salt bridges leave some pockets
electrostatically incomplete). Neither SQM nor AEV-PLIG can score the three silicon-rhodamine dye
systems (COSMO2 has no Si parameters; AEV-PLIG rejects the Si poses), so those cells are blank; GNINA
scores all 60.

## Repository layout

```
reference/system_reference.csv                      id, protein, ligand_name, sequence, smiles, source, source_url
reference/experimental_reference_ground_truth.csv   id, is_binder, experimental_pKD
predictions/<method>_predictions.csv                method, id, predicted_affinity
scripts/score.py, scripts/make_figure.py            the statistics, and the figure
scripts/methods.py                                  registry: every column, its producer, its category
scripts/append_paper.py                             the only thing that writes reference/
```

Every `predictions/*.csv` joins to the reference by `id` and is scored automatically. Adding a
column means adding a file and registering it; `scripts/qc_dataset.py` runs 21 integrity checks on
every PR.

## Extending the set

The set grows by a few papers a year. A weekly job searches Europe PMC, drafts an extraction and
opens a pull request; a person reviews and merges it, because a wrong row is worse than a missing
one. `scripts/append_paper.py` is the only thing that writes `reference/`.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the ingest loop, how to score a new row (hosted on
Rowan, or on a GPU for the 11 columns that need one), the workflows and what bounds their spending.

## Reproduce

```
# each cofolder needs its own GPU env; baselines are RDKit
CUDA_VISIBLE_DEVICES=0 BOLTZ=/path/to/boltz            python scripts/run_predictors.py boltz
CUDA_VISIBLE_DEVICES=0 NESSO=/path/to/nesso            python scripts/run_predictors.py nesso
CUDA_VISIBLE_DEVICES=0 PROTENIX_DIR=/path/to/Protenix  python scripts/run_predictors.py protenix
python scripts/run_predictors.py baselines
python scripts/run_predictors.py protenix-parse  # re-score a saved Protenix fold, no GPU
python scripts/run_ipdae.py --report   # ipDAE from the saved Boltz-2 PAE; needs _work/boltz_out
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
[Rowan Scientific](https://rowansci.com) ([@RowanSci](https://x.com/RowanSci)), with thanks to
[Corin Wagen](https://x.com/cwagen) for the credits. See [References](#references) for the models,
scorers, and libraries this benchmark builds on.
