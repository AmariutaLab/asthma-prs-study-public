# `sample_data/` — demo fixture for PTRS scoring (Stage 1)

These files let anyone run [`../run_demo.sh`](../run_demo.sh) end-to-end (the
actual `ptrs_score_groupC_v4.R` worker, tissue track) without the full GTEx
weight panel or the controlled-access study cohorts. **No individual-level study
data is included.**

Scope: **1 tissue (Whole_Blood) · chromosome 1 · 12 genes** — the same genes as
the [`05_twas-analysis/`](../../05_twas-analysis/) demo.

| File | Real or fake? | Provenance |
|------|---------------|------------|
| `weights/META_Whole_Blood/*.wgt.RDat` (12) | **REAL public data** | GTEx v8 EUR FUSION expression-weight models, Whole_Blood, 12 chr1 genes. The worker picks each gene's best model (top1/lasso/enet/blup) by `cv.performance`. |
| `TranscriptsInWhole_BloodModel.demo.txt` | derived | the 12 gene ids (versioned), row-aligned with the TWAS z column. |
| `Marginal_alphas_NEW_TWAS_Whole_Blood.demo.txt.gz` | **REAL (derived)** | per-gene TWAS z-scores from the [`05_twas-analysis/`](../../05_twas-analysis/) demo run (used to z-weight the summed score). |
| `focus_geneList_39/chr1_Whole_Blood_gene_list_p1.txt` | derived | selected-gene list consumed by the worker (version-stripped ENSG ids). |
| `expected_outputs/chr1_Whole_Blood_all_data{,_keep,_sign}.rds` | reference run | committed per-(chr, tissue) PTRS scores for regression comparison. |

The **genotype target** is not committed here — the runner reuses the public
1000 Genomes EUR chr1 fileset from
`../../05_twas-analysis/sample_data/LDREF/1000G.EUR.1` as a stand-in for the
controlled-access cohort target (GACRS/CAMP). The worker rewrites the target
`.bim` in place during allele harmonization, so the runner copies it into
`demo_run/` first and never touches the committed copy.

## Notes on faithfulness

- The runner invokes the **unmodified** `ptrs_score_groupC_v4.R`; the only
  scaffolding is a `plink` symlink under `PROJECT_DIR` (the worker calls
  `${PROJECT_DIR}/./plink`) and the target copy.
- The demo's gene list is the full 12-gene demo set (for a fuller illustration
  of the cross-gene sum). In production this list is the **FOCUS credible set**
  from [`07_focus-finemapping/`](../../07_focus-finemapping/) (`focus_geneList_39/`)
  or the [`06_TWAS_PT/`](../../06_TWAS_PT/) clumping output.
- Stage 2 (cross-chromosome concat + 80/20 cohort split) is the
  `PTRS_allChr_score_*_unified.RMD` notebooks; they consume these per-chromosome
  RDS files. Not run here (needs cohort phenotypes), but documented in the
  module README.

## Cell-type (OneK1K) track

Not included — the `v5` cell-type worker needs the OneK1K FUSION weight panels,
which are not available to redistribute here. See the module README,
*Cell-type (OneK1K) track — pending*.
