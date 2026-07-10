# 08 — PTRS Construction

End-to-end PTRS construction: from per-gene predicted expression to per-cohort
PTRS result CSVs consumed by the unified notebooks in
[`09_ptrs-unified_model-evaluation/`](../09_ptrs-unified_model-evaluation/).

> [!TIP]
> **Reproducible demo (smoke test).** A committed fixture + runner exercise the
> **unmodified** Stage-1 worker end-to-end on public data:
> ```bash
> cd 08_ptrs-construction
> PLINK=plink RSCRIPT=Rscript ./run_demo.sh
> ```
> It scores **1 GTEx tissue (Whole_Blood), chromosome 1, 12 genes** using real
> GTEx v8 FUSION weights + the TWAS z from the 05 demo + a public 1000G-EUR chr1
> genotype target — all in [`sample_data/`](sample_data/) (~110 KB). It writes
> the three per-(chr, tissue) RDS scores (raw / z-weighted / sign-weighted) and
> checks them against the committed reference. Requires R with `data.table`,
> `glmnet`, `dplyr` and PLINK 1.9 (v1.9.0-b.7.7). See
> [`sample_data/README.md`](sample_data/README.md).
>
> **Cell-type (OneK1K) track — pending.** The demo covers the tissue-track `v4`
> worker only; the `v5` cell-type worker needs the OneK1K FUSION weight panels,
> which are not redistributable here. **Stage 2** (the `*_unified.RMD`
> cross-chromosome concat + 80/20 cohort split) needs cohort phenotypes and is
> not part of the demo.

The pipeline has two stages:

1. **Per-(chromosome, tissue) scoring** — SLURM array job that imputes
   predicted expression for each significant gene (`plink --score`), weights it
   by the per-gene TWAS z-score, and sums across genes into three RDS files
   per chromosome.
2. **Cross-chromosome concatenation + selection** — RMDs that read those RDS
   files, sum `SCORE` across chromosomes 1-22, split the GACRS cohort
   (80/20, seed=42), and emit per-feature CSVs plus a validation-set focus
   list.

Stage 1 has two parallel script pairs, one per feature set:

- **v4** (`ptrs_score_v4.sbatch` + `ptrs_score_groupC_v4.R`) — 39 GTEx
  bulk-tissue variant. Reads the 6-column tissue TWAS path table (with a
  real `transcript_keep_file` filter), pulls gene lists from
  `focus_geneList_39/`, and supports `gacrs_1kg` / `camp_1k1k` cohorts.
- **v5** (`ptrs_score_v5.sbatch` + `ptrs_score_groupC_v5.R`) — 17 OneK1K
  cell-type variant. Reads the 5-column CT TWAS path table (no separate
  filter), pulls gene lists from `focus_geneList_17CT/`, and supports
  `gacrs_1kg` / `camp_gtex` cohorts.

The two `.R` workers are functionally identical after cleanup — both are
fully parameterized via 10 positional args — and are kept as separate files
purely for naming parity with the v4/v5 sbatch drivers.

## Environment setup

Stage 1 is **R** + PLINK; Stage 2 is the `*_unified.RMD` notebooks (R Markdown).

```bash
# R worker packages (the study used R 4.3.1)
Rscript -e 'install.packages(c("data.table", "glmnet", "dplyr"))'
# Stage 2 RMDs additionally use: rmarkdown, knitr

# per-gene scoring
conda create -n asthma-demo -c bioconda plink=1.90b7.7      # PLINK 1.9 v1.9.0-b.7.7
```

The worker calls `${PROJECT_DIR}/./plink`; `run_demo.sh` sets that up via a
symlink so the unmodified worker runs against a `plink` on your PATH.

## Tissue vs CT: heritability keep-list handling

Both tracks score only the heritability-QC-passing subset of genes; the two
just hit that filter at different points in the pipeline. See
`05_twas-analysis/README.md` — *Tissue vs CT: heritability keep-list
handling* — for the canonical side-by-side.

- **Tissue (`ptrs_score_groupC_v4*.R`)** reads the `transcript_keep_file`
  column from the 6-column tissue TWAS path table — it points at the
  pre-given TCSC keep file `TranscriptsIn<Tissue>Model_keep.txt`. The
  worker `merge`s the transcripts file against the keep file with
  `all.y = TRUE`, dropping genes outside the keep list before scoring.
  Tissue is the one stage where the keep list is re-applied at score time.
- **CT (`ptrs_score_groupC_v5*.R`)** reads from the 5-column CT TWAS path
  table, which has no `transcript_keep_file` column. None is needed: the
  `transcripts` column already points at a per-CT transcripts file that
  was intersected with the inline-generated `twas_qc_genelist_ct.csv` by
  `05_twas-analysis/Rebuild_twas_path_table_ct.R`. The CT worker therefore
  effectively scores the same kind of QC-passing subset that the tissue
  worker does — it just doesn't have to re-do the intersection here.

The same parallel holds for the TWAS P+T variant (`*_twas_pt`): the tissue
worker re-applies `_keep`, the CT worker doesn't, and both end up scoring
heritability-QC-passing genes only.

## Files

| File | Stage | Purpose |
|---|---|---|
| `ptrs_score_v4.sbatch` | 1 | Tissue (39 GTEx) SLURM array driver. Sets `COHORT` at the top, looks up per-tissue weight + transcript paths from the 6-column tissue TWAS path table, and invokes `ptrs_score_groupC_v4.R` once per tissue per chromosome. |
| `ptrs_score_groupC_v4.R` | 1 | Tissue-mode worker. Picks the best FUSION model (top1 / lasso / enet / blup) by `cv.performance` p-value, with the groupC single-column path. Strand-flips / allele-swaps SNPs to match the weight file's effect allele, runs `plink --score` per gene, and writes `chr<chr>_<tissue>_all_data{,_keep,_sign}.rds`. |
| `ptrs_score_v5.sbatch` | 1 | CT (17 OneK1K) SLURM array driver. Same shape as v4 but reads the 5-column CT TWAS path table and calls `ptrs_score_groupC_v5.R`. |
| `ptrs_score_groupC_v5.R` | 1 | CT-mode worker. Functionally identical to v4.R; see header. |
| `PTRS_allChr_score_gacrs_unified.RMD` | 2 | Concats per-chromosome RDS files for each feature, merges with GACRS phenotype, splits GACRS into 80% train+test / 20% validation (seed=42), writes per-feature CSVs, and runs validation-set feature selection to produce the focus list (`selected_*_all_case_gt_ctrl_keep.csv`). |
| `PTRS_allChr_score_camp_unified.RMD` | 2 | Same concat step on the external cohort: CAMP+1KG (tissue) or CAMP+GTEx (CT). No train/test split — entire cohort is external validation. |

## How to run

### Stage 1 — per-(chr, tissue) scoring (SLURM)

For tissue (39 GTEx):

```bash
# COHORT ∈ {gacrs_1kg, camp_1k1k}; defaults to gacrs_1kg
sbatch --export=ALL,COHORT=gacrs_1kg,PVAL=1 ptrs_score_v4.sbatch
```

For CT (17 OneK1K):

```bash
# COHORT ∈ {gacrs_1kg, camp_gtex}; defaults to camp_gtex
sbatch --export=ALL,COHORT=camp_gtex,PVAL=1 ptrs_score_v5.sbatch
```

Both jobs are 22-element arrays (one chromosome per task) and loop over every
tissue / cell type in the chosen feature set. Per-tissue paths (weights,
transcripts, TWAS z-scores) are looked up from the TWAS path table
(`twas_path_table_39_new.csv` for tissue, `twas_path_table_ct_new.csv` for ct).

The R workers accept 10 positional args; both `plink_basename` and
`gene_list_template` are derived in the sbatch from `COHORT`, so the R script
never needs to be edited to switch cohorts. The gene-list template supports
`{chr}`, `{tissue}`, `{pval}` placeholders, e.g.
`focus_geneList_17CT/fdr_sig_genes_chr{chr}_{tissue}.txt`.

### Stage 2 — concat + selection (R Markdown)

Each RMD is parameterized by a single `MODEL_VERSION` constant at the top
(`"tissue"` or `"ct"`). Flip the value in the first chunk and re-knit to
produce the other variant.

```r
# In RStudio: open the .RMD, set MODEL_VERSION in chunk 1, then Knit.
# From the shell:
Rscript -e 'rmarkdown::render("PTRS_allChr_score_gacrs_unified.RMD")'
Rscript -e 'rmarkdown::render("PTRS_allChr_score_camp_unified.RMD")'
```

To regenerate both variants you must run each file twice — once with
`MODEL_VERSION = "tissue"` and once with `MODEL_VERSION = "ct"`.

## Configuration matrix

| sbatch | `COHORT` | PLINK basename | Stage-1 OUTPUT_DIR | TWAS table | Gene-list dir |
|---|---|---|---|---|---|
| `ptrs_score_v4.sbatch` (tissue, 39) | `gacrs_1kg` | `GACRS_1kg/GACRS_1kg_chr` | `outputs_39_focus_gacrs_1kg_new` | `TWAS/rebuilt/twas_path_table_39_new.csv` (6 col) | `focus_geneList_39/` |
| `ptrs_score_v4.sbatch` (tissue, 39) | `camp_1k1k` | `CAMP_onek1k/CAMP_onek1k_chr` | `outputs_39_focus_camp_1k1k` | (same) | (same) |
| `ptrs_score_v5.sbatch` (CT, 17)     | `gacrs_1kg` | `GACRS_1kg/GACRS_1kg_chr` | `outputs_17CT_focus_gacrs_1kg` | `TWAS/rebuilt_ct/twas_path_table_ct_new.csv` (5 col) | `focus_geneList_17CT/` |
| `ptrs_score_v5.sbatch` (CT, 17)     | `camp_gtex` | `CAMP_GTEX/CAMP_GTEX_chr` | `outputs_17CT_focus_camp_gtex` | (same) | (same) |

## Outputs

### Stage 1 (per chromosome × tissue)

`<OUTPUT_DIR>/<pval>/<tissue>/chr<chr>_<tissue>_all_data{,_keep,_sign}.rds`

Each RDS is a data.table with `FID, IID, PHENO, CNT, CNT2, SCORE` where
`SCORE` is summed over the chromosome's significant genes:

- `all_data`      — raw sum of `plink --score` per-gene scores
- `all_data_keep` — weighted by signed gene TWAS z-score
- `all_data_sign` — weighted by `sign(z)` only

### Stage 2 (cross-chromosome concat + selection)

| Output | tissue (39) | ct (17) |
|---|---|---|
| Per-feature concat CSVs (GACRS train+test) | `combine/data/1/ptrs_results-concat_39_gacrs_train_test/<tissue>_results.csv` | `combine/data/1/ptrs_results-concat_17CT_gacrs_train_test/<ct>_results.csv` |
| Per-feature concat CSVs (GACRS val) | `combine/data/1/ptrs_results-concat_39_gacrs_val/...` | `combine/data/1/ptrs_results-concat_17CT_gacrs_val/...` |
| Per-feature concat CSVs (external) | `combine/data/1/ptrs_results-concat_39_camp_1k1k/...` | `combine/data/1/ptrs_results-concat_17CT_camp_gtex/...` |
| Validation-set selection (notebook input) | `results_concat_39_new/selected_gacrs_val_all_case_gt_ctrl_keep.csv` | `results_concat_17CT/selected_17CT_gacrs_val_all_case_gt_ctrl_keep.csv` |
| GACRS phenotype splits | `combine/data/updated_pheno_gacrs_{val,train_test}.txt` | (shared with tissue) |
| External cohort phenotype | `combine/data/updated_pheno_camp_1k1k.txt` | `combine/data/updated_pheno_camp_gtex.txt` |

## Inputs expected

**Stage 1 (both v4 and v5):**

- TWAS path table —
  - tissue: `TWAS/rebuilt/twas_path_table_39_new.csv` (6 columns:
    `tissue, gtex_twas, transcripts, transcript_keep_file, weights, PANEL`)
  - ct: `TWAS/rebuilt_ct/twas_path_table_ct_new.csv` (5 columns: same minus
    `transcript_keep_file`)
- Per-gene FUSION weight files: `<weight_file_basename><gene_id>.wgt.RDat`
- Per-tissue transcript list + matching TWAS z-score column
- Significant-gene lists per (chr, tissue) at the path produced by
  `gene_list_template`, e.g.
  `focus_geneList_17CT/fdr_sig_genes_chr<chr>_<tissue>.txt` or
  `focus_geneList_39/fdr_sig_genes_chr<chr>_<tissue>.txt`
- Cohort PLINK files (`<plink_basename><chr>.{bim,bed,fam}`).
  **Note:** Stage 1 rewrites the cohort BIM in place to align A1/A2 with the
  weight file's effect allele; keep a backup if you need the original.
- PLINK binary at `${PROJECT_DIR}/plink`

**Stage 2:**

- Stage-1 RDS files under
  `outputs_{39_focus,17CT_focus}_{gacrs_1kg,camp_1k1k,camp_gtex}/1/<feature>/chr<N>_<feature>_all_data{,_keep,_sign}.rds`
- GACRS: `GACRS_1kg_chr1.fam` + `gacrs_cases.txt` (case list)
- External: `CAMP_onek1k_chr1.fam` (tissue) or `CAMP_GTEX_chr1.fam` (ct)
- Tissue list (for `tissue` version): `twas_path_table_39.csv`

## Data inputs (real vs sample candidates)

These are the categories of input data this stage consumes. Sizes are
approximate; "Sample viable?" flags which inputs are small or anonymizable
enough that we could ship a committed example fixture vs which would have to
remain external. The cohort PLINK files contain identifiable genotype
data — these almost certainly cannot be shipped as-is. The final
real-vs-sample choice for each row is left for the maintainer to fill in.

| Input | Source | Approx. size | Sample viable? | Notes |
|---|---|---|---|---|
| TWAS path table (tissue, 6-col) | `05_twas-analysis/` rebuilt/ | ~KB | **Yes** | Committable |
| TWAS path table (CT, 5-col) | `05_twas-analysis/` rebuilt_ct/ | ~KB | **Yes** | Committable |
| Per-gene FUSION weight files (`<weight_basename><gene_id>.wgt.RDat`) | external (GTEx v8 / OneK1K) | ~10 GB total | **No** | Shared external dep with 05 / 07 |
| Per-tissue transcripts + TWAS Z-scores | `05_twas-analysis/` rebuilt{,_ct}/ | ~MB / feature | **Yes** | Same fixtures used by 06 |
| Per-(chr, tissue) gene lists from FOCUS (`focus_geneList_{39,17CT}/`) | `07_focus-finemapping/` | ~KB / file | **Yes** | A handful is committable; ship enough to cover one tissue × one chr |
| Per-(chr, tissue) gene lists from P+T (`gene_pt_{new,ct_new}/twas_gene_lists/`) | `06_TWAS_PT/` | ~KB / file | **Yes** | Same — ship a handful per p-value |
| Cohort PLINK files — GACRS + 1KG (`GACRS_1kg_chr<N>.{bim,bed,fam}`) | external (cohort) | ~5 GB / cohort | **No** | Identifiable; document data-access path. **Sub-option**: ship a small synthetic PLINK fixture (N=10) for end-to-end smoke test |
| Cohort PLINK files — CAMP+OneK1K, CAMP+GTEx | external (cohort) | ~5 GB / cohort | **No** | Same |
| GACRS phenotype + case list (`gacrs_cases.txt`) | external (cohort) | ~KB | **Yes if anonymized** | Or synthetic phenotype matched to the synthetic PLINK fixture |
| PLINK binary | external | install | N/A | Link to upstream install instructions |
| Stage-1 RDS outputs (for Stage 2 testing) | from Stage 1 | ~MB / (chr, feature) | **Yes** | A handful is committable as expected-output fixtures |

## TWAS P+T variant (`*_twas_pt`)

Stage 1 ships a second script pair that scores the **same** per-gene model but
draws its significant-gene lists from the TWAS p-value-thresholding (P+T)
pipeline — LD-clumped, per-pval gene lists — instead of the FOCUS fine-mapped
lists used by the `ptrs_score_v4/v5` scripts above. The FOCUS scripts are
unchanged; the P+T variant is purely additive.

| File | Mirrors | Difference |
|---|---|---|
| `ptrs_score_v4_twas_pt.sbatch` | `ptrs_score_v4.sbatch` | Tissue driver; `GENE_LIST_TEMPLATE` → `TWAS/gene_pt_new/twas_gene_lists/`, `OUTPUT_DIR` → `outputs_39_pt_*` |
| `ptrs_score_groupC_v4_twas_pt.R` | `ptrs_score_groupC_v4.R` | Tissue worker; identical scoring loop |
| `ptrs_score_v5_twas_pt.sbatch` | `ptrs_score_v5.sbatch` | CT driver; `GENE_LIST_TEMPLATE` → `TWAS/gene_pt_ct_new/twas_gene_lists/`, `OUTPUT_DIR` → `outputs_17CT_pt_*` |
| `ptrs_score_groupC_v5_twas_pt.R` | `ptrs_score_groupC_v5.R` | CT worker; resolves the gene-list `{tissue}` token with the case-preserved cell-type name (weights-dir basename, e.g. `B_intermediate`) — P+T CT gene-list filenames are cased |

The scoring loop (model selection, strand-flip / allele-swap, `plink --score`,
the three `all_data{,_keep,_sign}` accumulators) is identical to the FOCUS
workers — only the gene-list source and the output directory differ. The P+T
variant writes to its own `outputs_*_pt_*` directories, so P+T and FOCUS
results never overwrite each other.

### How to run

```bash
# tissue: COHORT ∈ {gacrs_1kg, camp_1k1k}; PVAL ∈ {5e-05, 5e-04, 0_005, 0_05}
sbatch --export=ALL,COHORT=gacrs_1kg,PVAL=0_005 ptrs_score_v4_twas_pt.sbatch

# ct:     COHORT ∈ {gacrs_1kg, camp_gtex};  PVAL ∈ {5e-05, 5e-04, 0_005, 0_05}
sbatch --export=ALL,COHORT=gacrs_1kg,PVAL=0_005 ptrs_score_v5_twas_pt.sbatch
```

**Concurrency note:** the tissue and CT GACRS jobs both rewrite the same
`GACRS_1kg/*.bim` in place. Do not run them at the same time — chain them with
`--dependency=afterany` so only one job touches the shared BIM at a time.

### Configuration matrix (P+T)

| sbatch | `COHORT` | OUTPUT_DIR | Gene-list dir |
|---|---|---|---|
| `ptrs_score_v4_twas_pt.sbatch` (tissue, 39) | `gacrs_1kg` | `outputs_39_pt_gacrs_1kg_new` | `TWAS/gene_pt_new/twas_gene_lists/` |
| `ptrs_score_v4_twas_pt.sbatch` (tissue, 39) | `camp_1k1k` | `outputs_39_pt_camp_1k1k_new` | (same) |
| `ptrs_score_v5_twas_pt.sbatch` (CT, 17)     | `gacrs_1kg` | `outputs_17CT_pt_gacrs_1kg_new` | `TWAS/gene_pt_ct_new/twas_gene_lists/` |
| `ptrs_score_v5_twas_pt.sbatch` (CT, 17)     | `camp_gtex` | `outputs_17CT_pt_camp_gtex_new` | (same) |

P+T gene-list filenames follow `chr<chr>_<tissue>_gene_list_p<pval>.txt` — CT
uses the case-preserved cell-type name, e.g.
`chr1_B_intermediate_gene_list_p0_005.txt`.
