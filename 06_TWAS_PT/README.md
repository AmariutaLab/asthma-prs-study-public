# 06_TWAS_PT — TWAS p-value LD-clumping gene selection

Alternative gene-selection pathway for PTRS construction. Where
`07_focus-finemapping` uses ma-FOCUS posteriors to pick a credible gene
set, this directory picks genes by **TWAS p-value rank with LD clumping**
on predicted expression — a simpler, faster baseline that the downstream
PTRS pipeline (`08_ptrs-construction`) consumes in the same shape.

Pipeline position: TWAS Z-scores (`05_twas-analysis`) → **06_TWAS_PT (this dir)** → PTRS scoring (`08_ptrs-construction`).

> [!TIP]
> **Reproducible demo (smoke test).** A committed fixture + runner exercise the
> tissue-track clumping method end-to-end on public data:
> ```bash
> cd 06_TWAS_PT
> RSCRIPT=Rscript ./run_demo.sh
> ```
> It runs the documented method (correlation → p from Z → sort → greedy
> r²<0.1 clump → keep-filter → threshold + per-chr split) on **Whole_Blood, 12
> chr1 genes**, using a real GTEx v8 predicted-expression subset + the TWAS Z
> from the 05 demo — all in [`sample_data/`](sample_data/) (~40 KB). Requires R
> with `data.table`. See [`sample_data/README.md`](sample_data/README.md).
>
> **Cell-type (OneK1K) track — pending.** The demo covers the tissue track only;
> the cell-type track needs the OneK1K predicted-expression CSVs, which are not
> redistributable here.

## Environment setup

Pure **R** (the study used R 4.3.1):

```r
install.packages(c("data.table", "glmnet"))
```

## Inputs

- TWAS Z-scores per tissue / cell type
  - Tissue: `Marginal_alphas_NEW_TWAS_<tissue>.txt.gz` (39 GTEx tissues)
  - Cell type: `Marginal_alphas_NEW_TWAS_<celltype>.txt.gz` (17 OneK1K CTs)
  - One column of Z-scores, row-aligned with the matching `TranscriptsIn<*>Model.txt`
- Predicted expression matrices
  - Tissue: GTEx v8 `designmat_<tissue>_v8_{320EUR,allEUR}_double.RData` (binary `df1`, columns labeled separately by `TranscriptsIn<tissue>Model.txt`)
  - Cell type: `<celltype>_pred_expr.csv` from collaborator (header carries gene symbols)
- Gene annotation: `gene_annotation.txt.gz` (chr, start, stop, symbol, …, Ensembl)

## Tissue vs CT: heritability keep-list handling

Both tracks apply the same per-feature heritability whitelist; the two
just differ in whether the whitelist is re-applied here or has already been
baked into the upstream input. See `05_twas-analysis/README.md` —
*Tissue vs CT: heritability keep-list handling* — for the canonical
side-by-side.

- **Tissue (`ld_clumping_tissue.R`)** clumps by predicted-expression
  r² < 0.1 and then intersects with the pre-given TCSC keep list
  `TranscriptsIn<tissue>Model_keep.txt` as a post-clump step.
- **CT (`ld_clumping_celltype.R`)** clumps by predicted-expression
  r² < 0.1 and stops there. No intersection with `twas_qc_genelist_ct.csv`
  is needed because the inputs the CT clumping reads (the TWAS Z file and
  the matching predicted-expression CSV) have already been restricted to
  QC-passing genes upstream in `05_twas-analysis/Rebuild_twas_path_table_ct.R`.

Functionally the two pipelines apply the same heritability whitelist —
tissue does it twice (once upstream when building the path table, again
here as the post-clump intersection), and CT does it once (upstream only).

## Method

1. Build a gene-gene correlation matrix from each tissue/CT's predicted-expression matrix (`cor()`).
2. Convert TWAS Z to a two-sided p-value (`p = 2 * pnorm(|Z|, lower.tail = FALSE)`).
3. Sort genes by p ascending.
4. Greedy LD clump: walk in p-order; keep a gene iff r² < 0.1 against every already-selected gene.
5. (Tissue only) Intersect with the heritability keep-list `TranscriptsIn<tissue>Model_keep.txt`. CT pipeline has no keep-list step — QC is enforced upstream.
6. Emit gene lists at four p-thresholds: `5e-5, 5e-4, 5e-3, 0.05`.
7. Split each list into per-chromosome files for the downstream array-job PTRS scoring.

## Files

| File | Purpose |
|---|---|
| `ld_clumping_tissue.R`   | Tissue clumping (39 GTEx tissues, Ensembl IDs, designmat input). |
| `ld_clumping_celltype.R` | Cell-type clumping (17 OneK1K CTs, HGNC symbols, collaborator CSV input). |
| `split_tissue.R`         | Fan each tissue gene list into per-chromosome files. |
| `split_celltype.R`       | Fan each cell-type gene list into per-chromosome files. |

## Outputs

`twas_gene_lists/` in the run directory:

- `<context>_gene_list_p<pval>.txt`
  - Row schema: `gene_id CHR START STOP <annotation>` (space-separated)
- `chr<N>_<context>_gene_list_p<pval>.txt`
  - One gene per line (the downstream ID — Ensembl base for tissue, symbol for CT)
  - Consumed by `08_ptrs-construction/ptrs_score_groupC_v{4,5}*.R`

`<context>` = tissue name (e.g., `Whole_Blood`) or case-preserved CT name
(e.g., `CD14_Mono`). `<pval>` = `5e-05`, `5e-04`, `0_005`, `0_05` (dots
replaced with underscores to keep shell-safe filenames).

## Run

```r
# tissue
Rscript ld_clumping_tissue.R
Rscript split_tissue.R

# cell type
Rscript ld_clumping_celltype.R
Rscript split_celltype.R
```

Each clumping script writes a `ld_clumping.log` next to its output dir.

## Data inputs (real vs sample candidates)

These are the categories of input data this stage consumes. Sizes are
approximate; "Sample viable?" flags which inputs are small enough to ship a
committed example fixture vs which would have to remain external downloads.
The final real-vs-sample choice for each row is left for the maintainer to
fill in.

| Input | Source | Approx. size | Sample viable? | Notes |
|---|---|---|---|---|
| Tissue TWAS Z-scores (`Marginal_alphas_NEW_TWAS_<tissue>.txt.gz`) | `05_twas-analysis/` rebuilt/ | ~MB / tissue | **Yes** | One tissue is committable; whole set is not too large either |
| Tissue transcripts (`TranscriptsIn<tissue>Model.txt`) | `05_twas-analysis/` rebuilt/ | ~KB / tissue | **Yes** | Tiny; commit alongside the matching Z-score file |
| CT TWAS Z-scores (`Marginal_alphas_NEW_TWAS_<ct>.txt.gz`) | `05_twas-analysis/` rebuilt_ct/ | ~MB / CT | **Yes** | One CT is committable |
| CT transcripts (`TranscriptsIn<ct>Model.txt`) | `05_twas-analysis/` rebuilt_ct/ | ~KB / CT | **Yes** | Commit alongside |
| GTEx predicted-expression matrices (`designmat_<tissue>_v8_*_double.RData`) | external (GTEx v8 / TCSC release) | ~50 MB / tissue | **Maybe** | One tissue is committable; whole set is ~2 GB |
| Cell-type predicted-expression CSVs (`<celltype>_pred_expr.csv`) | external collaborator | ~10 MB / CT | **Yes** | One CT is committable |
| Tissue heritability keep-list (`TranscriptsIn<tissue>Model_keep.txt`) | external (TCSC) | ~50 MB total | **Yes** | Same source as 05_twas-analysis |
| Gene annotation (`gene_annotation.txt.gz`) | shared with 05 | ~5 MB | **Yes** | Committable |

## Notes

- **Two-sided p**: cutoffs (`5e-5, 5e-4, 5e-3, 0.05`) are applied to `2 * pnorm(|Z|, lower.tail = FALSE)`.
- **Clumping is order-dependent**: rank-equivalent to |Z|, so sort order is invariant to one-vs-two-sided p; only the threshold cutoff changes.
- **No transcripts-keep filter for CT**: the cell-type path table has no `_keep.txt` column. Tissue retains the heritability filter as a post-clump intersection.
- Gene-name conventions: tissue clumping strips Ensembl version (`ENSG…\.\d+` → `ENSG…`); CT uses raw HGNC symbols (no version).
