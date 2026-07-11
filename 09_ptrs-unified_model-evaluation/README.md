# 09 — Unified PRS + PTRS Model Evaluation

Per-feature PTRS exploration (MA-FOCUS and TWAS P+T), PRS-CS/PRS-CSx
evaluation, and the cross-modal integrated PRS + PTRS model evaluated on
CAMP-only Balanced (64v64 × 100 bootstrap). All notebooks are hard-coded to
the **altPRS** config (PRS-CS ϕ=auto/EUR, PRS-CSx ϕ=auto/META) and the
CAMP-only cohort. External-cohort (`CAMP+GTEx`, `CAMP+1KG`) columns and
Nagelkerke R² are no longer reported.

> [!TIP]
> **Reproducible demo (smoke test).** The notebooks run on individual-level
> per-sample predictions (`data/`), which are not public. A compact,
> **fully-synthetic** stand-in reproduces the headline computation end-to-end:
> ```bash
> cd 09_ptrs-unified_model-evaluation
> python run_demo.py            # needs scikit-learn, pandas, numpy (the `test` env)
> ```
> [`run_demo.py`](run_demo.py) mirrors `integrated_ptrs_prs_unified.ipynb` —
> per-feature **OOF transform** → **Direct integration** of
> `cd4_naive + Esophagus_Mucosa` PTRS with PRS → **CAMP-Balanced bootstrap AUC**
> vs PRS-alone — on the synthetic fixture in [`sample_data/`](sample_data/). It
> reproduces the paper's *ordering and mechanism* (integrated > PTRS-only >
> PRS-alone; ΔAUC > 0), is deterministic (matches the committed
> `expected_outputs/` exactly), and uses **no individual-level data** — the real
> `data/` directory is untouched. See [`sample_data/README.md`](sample_data/README.md).

## Environment setup

Python (scikit-learn stack); the study ran these in the `test` conda env.

```bash
conda create -n ptrs-eval -c conda-forge \
    python=3.11 scikit-learn pandas numpy matplotlib jupyter
conda activate ptrs-eval
```

`run_demo.py` needs only `scikit-learn`, `pandas`, `numpy`; the notebooks add
`matplotlib`/`jupyter` for the figure and summary steps.

## Pipeline overview

Four notebooks chain together to build and evaluate the integrated model:

1. **`meta_model_exploration_unified.ipynb`** — MA-FOCUS per-feature meta-model
   exploration for both modalities (39 GTEx tissues / 17 OneK1K cell types).
   Selects each feature's best classifier and writes the cross-cohort
   `consistent_features.csv` shortlist consumed downstream.

2. **`meta_model_exploration_TWAS_PT.ipynb`** — parallel TWAS P+T per-feature
   exploration at four LD-clumping p-value thresholds
   (`5e-05`, `5e-04`, `0.005`, `0.05`). Writes per-p-value
   `meta_model_{tissue,ct}__pval-*/` result directories used by the
   TWAS-P+T rows in the summary table + Figure 4B.

3. **`prscs_evaluation.ipynb`** — PRS-CS and PRS-CSx evaluation. Hard-coded to
   the **altPRS** configuration that was validated as best for the integrated
   model:
   - **PRS-CS** at ϕ=auto (LDREF=EUR, single-population)
   - **PRS-CSx** at ϕ=auto, LDREF=META (cross-population calibration)

4. **`integrated_ptrs_prs_unified.ipynb`** — the rank-1 integrated model:
   **cross-modal** `cd4_naive + Esophagus_Mucosa` PTRS with **per-feature OOF
   transform**, integrated directly with the altPRS PRS columns via 6
   classifiers + Rank Addition (7 total, matching the Methods "Integration
   strategies" section — Stacking-Tuned is intentionally omitted).

A fifth notebook, **`summary_table_and_plots.ipynb`**, reads the saved
predictions from all four and produces the cross-category comparison table
(`summary_comparison_table.csv`), consumed by `build_supplementary_tables.py`.
(Its former violin + pairwise-P outputs were removed — those figures are now
built by the Figure 4 pipeline in [`figure_scripts/`](figure_scripts/).)

## Notebooks

| File | Purpose |
|---|---|
| `meta_model_exploration_unified.ipynb` | **MA-FOCUS** per-feature PTRS + 10-method unified-PTRS exploration; writes `consistent_features.csv` shortlist per modality |
| `meta_model_exploration_TWAS_PT.ipynb` | **TWAS P+T** per-feature + unified exploration at four LD-clumping p-value thresholds (`5e-05`, `5e-04`, `0.005`, `0.05`); writes `meta_model_{tissue,ct}__pval-*/` directories used by the TWAS-P+T rows in the summary table |
| `prscs_evaluation.ipynb` | PRS-CS / PRS-CSx evaluation with train / val / test split; hard-coded to altPRS configs (PRS-CS ϕ=auto, PRS-CSx ϕ=auto/META); saves per-sample predictions to `data/predictions/prscs_evaluation/prs_predictions.csv` |
| `integrated_ptrs_prs_unified.ipynb` | **Cross-modal `cd4_naive + Esophagus_Mucosa` integrated with PRS via per-feature OOF transform + Direct integration**; 6 model-based classifiers (LR · Elastic Net CV · RF tuned · GB tuned · SVM linear · Stacking-Fixed LR+RF+GB) + Rank Addition = 7 total, matching Methods "Integration strategies". Stacking-Tuned is intentionally omitted. Evaluated on **CAMP-only Balanced 64v64 × 100 bootstrap** (no external cohort). |
| `integrated_ptrs_prs_unified_unimodal.ipynb` | **Uni-modal unified PTRS + PRS** integration. `MODEL_VERSION ∈ {'tissue', 'ct'}` switch in cell 0 selects the modality; runs the same 7-classifier zoo + Direct (PRS-only) baseline on the matched `consistent_features` shortlist for that modality, integrated with PRS-CS / PRS-CSx (altPRS). Writes per-modality results under `data/predictions/integrated_ptrs_prs_{tissue,ct}/`. Rank-1 uni-modal rows now use altPRS + 7-classifier grid to match cross-modal. |
| `summary_table_and_plots.ipynb` | Cross-category comparison table (`summary_comparison_table.csv`), ranked by CAMP-only Balanced AUC; consumed by `build_supplementary_tables.py`. (Obsolete violin + pairwise-P plotting removed — see the Figure 4 pipeline.) |
| `feature_count_ablation.py` | **Feature-count ablation** for the cross-modal architecture. Tests whether adding a third consistent-shortlist feature to the 2-feature anchor pair (cd4_naive + Esophagus_Mucosa) + PRS improves CAMP-Balanced AUC. Candidate list is read dynamically from `consistent_features.csv` — currently 9 candidates (7 GTEx tissues + 2 OneK1K cell types), 18 three-feature experiments × 2 PRS variants. Uses the paper's full 108-combo RF grid by default (`--reduced-grid` for a fast sanity check). Writes results and log to `data/predictions/feature_ablation/`. Reproduces the manuscript's directional claim: adding a 3rd feature drops CAMP-Balanced AUC in **18/18** experiments (mean ΔAUC = −0.031). |
| `compute_pairwise_bootstrap_pvalue.py` | **Pairwise bootstrap ΔAUC + one-sided P** helper. Reconstructs the 100-iteration CAMP-Balanced AUC series for each of 12 Figure-4B/4C models using the same `np.random.RandomState(i)` case-subsampling the notebooks use, then emits every ordered pair to `data/predictions/pairwise_comparisons/pairwise_bootstrap_auc.csv` (132 rows), plus the raw `per_iteration_auc_matrix.csv` (100 × 12) and per-model `model_mean_auc.csv`. Reproduces the paper's cited `P = 0.13` for Cross-modal + PRS-CSx RF vs Esophagus_Mucosa GB. |

## Headline result

The rank-1 model from `integrated_ptrs_prs_unified.ipynb`:

```
Direct (PRS-CSx) + Random Forest (tuned)
  Inputs:  [cd4_naive_OOF_prob, Esophagus_Mucosa_OOF_prob, PRS_CSx_z]
  CAMP-bal AUC = 0.6320 ± 0.040  (64 cases × 64 controls × 100 bootstrap)
  OR           = 3.55
```

The Figure-4C pairwise bootstrap comparisons anchored on this rank-1 model
(computed by `figure_scripts/build_figure4_panels.py`,
Panel C):

| Comparison | ΔAUC | one-sided P |
|---|---|---|
| Rank-1 vs PRS-CS alone | +0.141 | 0.00 |
| Rank-1 vs PRS-CSx alone | +0.119 | 0.03 |
| Rank-1 vs FOCUS-CT cd4_naive GB | +0.053 | 0.07 |
| Rank-1 vs FOCUS-tissue Esophagus_Mucosa GB | **+0.036** | **0.13** ← Fig 4C caption |
| Rank-1 vs Cross-modal + PRS-CS RF | +0.005 | 0.33 |

## Why cross-modal + per-feature OOF

- **Within-modality alone caps at ~0.59 on CAMP-bal**. The architecture-search
  in this folder tested every consistent feature individually, weighted sums,
  full unified-PTRS RF GridSearch, multiple PRS configs, and the 5-panel
  PRS-CSx ancestry breakdown — none crossed 0.60 within a single modality.
- **Cross-modal (1 cell-type + 1 tissue) breaks 0.60** because cd4_naive
  (immune-cell signal) and Esophagus_Mucosa (epithelial-tissue signal) are
  biologically distinct → genuinely orthogonal signal that combines
  constructively.
- **Per-feature OOF transform adds another +0.04 AUC** by calibrating each
  feature's raw `Keep_Vector` through its declared best per-feature model
  (RF GridSearch for cd4_naive, Gradient Boosting for Esophagus_Mucosa) — the
  unified-PTRS RF reduction step is **removed** because for 2 features it
  adds nothing over feeding both raw into the final classifier (Direct).
- **2 features is the sweet spot**. Adding any 3rd consistent feature to the
  anchor pair dropped CAMP-bal in **18/18** feature-count ablations
  (`feature_count_ablation.py`, 9 candidates × 2 PRS variants, full 108-combo
  RF grid). Mean ΔAUC = −0.031. More features → more dilution.

## Summary comparison table (`summary_table_and_plots.ipynb`)

Reads the saved prediction CSVs under `data/predictions/` and writes to `figures/`:

- `summary_comparison_table.csv` — best model per category (PRS-CS only,
  PRS-CSx only, best MA-FOCUS single-feature PTRS, best TWAS P+T
  single-feature PTRS, unified PTRS per modality, cross-modal integrated)
  with classifier, CAMP-only Balanced AUC, and ΔAUC vs the PRS-CSx baseline.
  Consumed by `build_supplementary_tables.py`.

> The notebook's former plotting half (`summary_violin1_prs_ptrs`,
> `summary_violin2_crossmodal`, and `summary_pairwise_vs_rank1.csv`) has been
> **removed** — those figures never entered the manuscript and are superseded by
> the Figure 4 pipeline (`build_figure4_panels.py` renders the violins as Panel B
> and the pairwise-P as Panel C).

Run from this folder — `ROOT = Path.cwd()` resolves all relative paths.

## Per-run outputs

`integrated_ptrs_prs_unified.ipynb` writes to
`data/predictions/integrated_ptrs_prs_combined/`:

- `per_feature_oof.csv` — per-sample OOF probabilities for each cross-modal
  feature (`cd4_naive`, `Esophagus_Mucosa`), across train OOF / test / CAMP
- `direct_predictions.csv` — per-sample integrated PRS + PTRS scores for every
  (classifier × PRS variant × eval set)
- `all_results.csv` — aggregated metric table (AUC, OR, P-value per row)

`integrated_ptrs_prs_unified_unimodal.ipynb` (run separately with
`MODEL_VERSION = 'tissue' | 'ct'`) writes to
`data/predictions/integrated_ptrs_prs_{tissue,ct}/`:

- `unified_ptrs.csv` — per-sample unified-PTRS score from the modality's
  `consistent_features` shortlist (RF-GridSearch reduction)
- `integrated_predictions.csv` — per-sample integrated PRS + unified-PTRS
  scores for every (classifier × PRS variant × eval set)
- `direct_predictions.csv` — per-sample PRS-only baseline scores from the
  same pipeline (used as the within-integration reference for ΔAUC)
- `all_results.csv` — aggregated metric table; rank-1 rows are
  RF (tuned) + PRS-CS for tissue (CAMP-bal AUC = 0.5302) and RF (tuned) +
  PRS-CS for CT (CAMP-bal AUC = 0.5031); rank-1 by AUC after the 7-classifier
  altPRS refactor.

`prscs_evaluation.ipynb` writes to `data/predictions/prscs_evaluation/`:

- `prs_predictions.csv` — per-sample PRS-CS / PRS-CSx scores at the altPRS
  config across GACRS Test / CAMP-only

`meta_model_exploration_unified.ipynb` (run separately per `MODEL_VERSION =
'tissue' | 'ct'`) writes to `data/predictions/meta_model_<MODEL_VERSION>/`:

- `consistent_features.csv` — feature shortlist consumed by the integrated
  notebook
- `individual_ptrs_long.csv`, `individual_3way_results.csv`, etc. — per-feature
  metric tables

## Input data (`data/`)

Per-cohort concat PTRS result CSVs used as inputs:

| Directory | Cohort split | Feature set |
|---|---|---|
| `ptrs_results-concat_39_gacrs_train_test/` | GACRS train+test | 39 GTEx tissues |
| `ptrs_results-concat_39_gacrs_val/`        | GACRS val        | 39 GTEx tissues |
| `ptrs_results-concat_39_train_new/`        | GACRS train      | 39 GTEx tissues |
| `ptrs_results-concat_39_val_new/`          | GACRS val        | 39 GTEx tissues |
| `ptrs_results-concat_39_camp_1k1k/`        | CAMP+1KG external | 39 GTEx tissues |
| `ptrs_results-concat_17CT_gacrs_train_test/` | GACRS train+test | 17 OneK1K cell types |
| `ptrs_results-concat_17CT_gacrs_val/`        | GACRS val        | 17 OneK1K cell types |
| `ptrs_results-concat_17CT_camp_gtex/`        | CAMP+GTEx external | 17 OneK1K cell types |

The integrated notebook only consumes the two `*_train_test/` directories +
the two external (`camp_gtex/` and `camp_1k1k/`) — it reads exactly one
per-feature CSV per directory (`cd4_naive_results.csv` and
`Esophagus_Mucosa_results.csv`).

Each subdirectory contains one CSV per feature
(`<tissue_or_cell_type>_results.csv`) with columns `Sample_ID`, `Keep_Vector`
(the PTRS score), and `asthma` (case/control label).

## Data inputs (real vs sample candidates)

These are the categories of input data the notebooks consume. Sizes are
approximate; "Sample viable?" flags which inputs are small or anonymizable
enough that we could ship a committed example fixture.

| Input | Source | Approx. size | Sample viable? | Notes |
|---|---|---|---|---|
| Per-cohort per-feature PTRS CSVs (`data/ptrs_results-concat_*/`) | `08_ptrs-construction/` Stage 2 | ~MB total (39 + 17 features × 3-4 cohorts) | **Yes** | Already lives in `data/`; subsample samples to N=100 to ship as fixtures, or keep real samples if anonymized. The integrated notebook only uses `cd4_naive_results.csv` from CT dirs and `Esophagus_Mucosa_results.csv` from tissue dirs. |
| PRS-CS / PRS-CSx melted CSVs (`${REPO_ROOT}/files/03_*_melt-admixture.csv`) | `04_prscs-prscsx-construction/` | ~10 MB total (6 files at ϕ=auto only) | **Yes** | Ship the 6 melt CSVs at ϕ=auto; the altPRS config only needs ϕ=auto rows |
| Phenotype files (`updated_pheno_gacrs_*`, `updated_pheno_camp_*`) | external (cohort metadata + 08 split) | ~KB | **Yes if anonymized** | Replace `IID` with synthetic IDs or commit a fake-sample fixture matched to the subsampled PTRS CSVs |
| FAM files for sample-to-cohort mapping (`GACRS_1kg_chr1.fam`, `CAMP_*_chr1.fam`) | external (cohort) | ~MB | **Yes if anonymized** | Same as above; can be regenerated from synthetic IDs |
| Saved fixtures (`data/predictions/meta_model_*/consistent_features.csv`) | from running `meta_model_exploration_unified.ipynb` | ~KB | **Yes** | Tiny CSVs the integrated notebook consumes as the feature-shortlist source |
| TWAS P+T fixtures (`data/predictions/meta_model_*__pval-*/`) | from running `meta_model_exploration_unified.ipynb` on the P+T inputs | ~KB / run | **Yes** | Already committed |
