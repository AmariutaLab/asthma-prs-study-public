# 05 — TWAS Analysis

FUSION TWAS over the meta-analyzed asthma GWAS, run twice in parallel:

- **39 GTEx bulk tissues** (tissue track), and
- **17 OneK1K cell types** (CT track).

> [!TIP]
> **Reproducible demo (smoke test).** A committed fixture + runner exercise the
> tissue-track FUSION association step end-to-end on public data:
> ```bash
> cd 05_twas-analysis
> FUSION_DIR=/path/to/fusion_twas RSCRIPT=Rscript ./run_demo.sh
> ```
> It runs `FUSION.assoc_test.R` on **1 GTEx tissue (Whole_Blood), chromosome 1,
> 12 genes**, using real GTEx v8 FUSION weights + a 1000G-EUR chr1 LD window +
> the asthma meta-GWAS sumstats sliced to those SNPs — all in
> [`sample_data/`](sample_data/) (~300 KB). The demo `TWAS.Z` matches the
> full-panel production run to <1e-2. Requires R with `plink2R`/`optparse` and a
> clone of [`gusevlab/fusion_twas`](https://github.com/gusevlab/fusion_twas). See
> [`sample_data/README.md`](sample_data/README.md) for provenance.
>
> **Cell-type (OneK1K) track — pending.** The demo covers the tissue track only.
> The 17-cell-type track needs the OneK1K eQTL weight panels, which are not
> redistributable here; that track's demo is deferred until those weights are
> obtained. The tissue-track fixture also stands in for the downstream demos in
> [`06_TWAS_PT/`](../06_TWAS_PT/), [`07_focus-finemapping/`](../07_focus-finemapping/),
> and [`08_ptrs-construction/`](../08_ptrs-construction/).

Outputs of this folder feed two downstream steps:

- The PTRS-construction stage 1 in
  [`08_ptrs-construction/`](../08_ptrs-construction/) consumes the rebuilt
  path tables (`rebuilt/twas_path_table_39_new.csv` and
  `rebuilt_ct/twas_path_table_ct_new.csv`) plus the row-aligned
  transcripts / TWAS-z files written next to them.
- The FOCUS fine-mapping step consumes
  `FOCUS/sorted_gene_all{39Tissue,17CT}_new.csv`.

## Environment setup

FUSION runs under **R** (the study used R 4.3.1); the LD reference + toolkit are
separate downloads.

```bash
# R packages for FUSION.assoc_test.R
Rscript -e 'install.packages(c("optparse", "RColorBrewer"))'
Rscript -e 'if(!requireNamespace("remotes",quietly=TRUE)) install.packages("remotes"); \
            remotes::install_github("gabraham/plink2R/plink2R")'   # reads LDREF/weights

# FUSION toolkit
git clone https://github.com/gusevlab/fusion_twas      # FUSION.assoc_test.R
```

The 1000G EUR `LDREF` panel and the GTEx/OneK1K FUSION weights are large public
downloads — see **Inputs expected**. (For the committed demo, `run_demo.sh` uses
a tiny bundled LD window + weights.)

## Pipeline

```
                  ┌──── TWAS_pos/      ─────────┐
                  │  (39 tissues, externally    │
[meta-GWAS] ──────┤   provided -- TODO)         ├──► [run FUSION assoc] ──► TWAS_results/
                  │                             │     submit_all_twas.sh
                  └─────────────────────────────┘             │
                                                              ├──► Rebuild_twas_path_table_39.R
                                                              │    └─► rebuilt/
                                                              │        ├─ gtex_twas/Marginal_alphas_NEW_TWAS_<T>.txt.gz
                                                              │        ├─ transcripts/TranscriptsIn<T>Model.txt
                                                              │        └─ twas_path_table_39_new.csv
                                                              │
                                                              └──► Make_sorted_gene_from_TWAS_results.R
                                                                   └─► FOCUS/sorted_gene_all39Tissue_new.csv

                  ┌──── TWAS_pos_ct/ ──────────┐
[GTEx celltype    │  Make_Pos_File_celltype.R  │
 weights v2]  ────┤                            ├──► [run FUSION assoc] ──► TWAS_results_ct/
                  │  (17 OneK1K cell types)    │     submit_all_twas_ct.sh
                  └─────────────────────────────┘             │
                                                              ├──► QC_twas_ct.R
                                                              │    ├─ twas_qc_genelist_ct.csv
                                                              │    └─ twas_qc_summary_ct.csv
                                                              │
                                                              ├──► Rebuild_twas_path_table_ct.R
                                                              │    └─► rebuilt_ct/
                                                              │        ├─ ct_twas/Marginal_alphas_NEW_TWAS_<ct>.txt.gz
                                                              │        ├─ transcripts/TranscriptsIn<ct>Model.txt
                                                              │        └─ twas_path_table_ct_new.csv
                                                              │
                                                              └──► Make_sorted_gene_from_TWAS_results_ct.R
                                                                   └─► FOCUS/sorted_gene_all17CT_new.csv
```

## Tissue vs CT: heritability keep-list handling

Both tracks apply the **same logical step** — restrict to genes whose
expression-prediction model passes heritability QC — using a per-feature gene
whitelist that is `merge`d / `intersect`ed in directly. The only difference
is **where that whitelist comes from**:

| | Tissue (39 GTEx) | CT (17 OneK1K) |
|---|---|---|
| Whitelist file | `TranscriptsIn<Tissue>Model_keep.txt` | `twas_qc_genelist_ct.csv` |
| Source | **Pre-given** — external TCSC heritable-genes lists at `${REF_DIR}/TCSC/weights/heritablegenes/{N320,Nall}/` | **Inline-generated** — `QC_twas_ct.R` writes it from `TWAS_results_ct/*.dat` using `HSQ > 0` + non-zero predictor variance + `MODELCV.R2 > 0` |
| Why the source differs | TCSC ships per-tissue heritable-genes lists with the GTEx panel, so we use them as-is | OneK1K has no equivalent pre-given list, so the same heritability filter is computed here from the TWAS `.dat` outputs |
| Where applied in 05 | `Make_Pos_File_tissue.R`, `Make_sorted_gene_from_TWAS_results.R`, and `Rebuild_twas_path_table_39.R` all intersect with the keep file | `Make_sorted_gene_from_TWAS_results_ct.R` and `Rebuild_twas_path_table_ct.R` both intersect with `twas_qc_genelist_ct.csv` |
| TWAS path table | 6 columns — explicit `transcript_keep_file` column points back at the per-tissue keep file, because the tissue worker in 08 re-applies it at score time | 5 columns — no `transcript_keep_file`, because the `transcripts` column has **already** been intersected with the QC genelist by `Rebuild_twas_path_table_ct.R` |
| Sorted-gene table for FOCUS (`sorted_gene_all*.csv`) | Built from the keep-filtered intersection | Built from the QC-genelist-filtered intersection |

Functionally, the two whitelists do the same job — they are not different
filters. The practical consequence for downstream stages (06_TWAS_PT,
07_focus-finemapping, 08_ptrs-construction) is just that **tissue carries
its keep file forward explicitly** (path-table column + re-applied at score
time in 08) while **CT bakes the equivalent filter into the artefacts it
emits** and no separate `_keep` reference is needed downstream.

## Files

| Stage | File | Track | Reads | Writes |
|---|---|---|---|---|
| 1. Build `.pos` | **(TODO)** | tissue | GTEx FUSION weights + gene annotation | `TWAS_pos/twas_<tissue>.pos` (externally provided; producer script will be added later) |
| 1. Build `.pos` | `Make_Pos_File_celltype.R` | CT | OneK1K celltype-weights v2 + gene annotation + per-CT `.fam` files | `TWAS_pos_ct/twas_<ct>.pos` |
| 2. Run FUSION TWAS | `submit_all_twas.sh` | tissue | `TWAS_pos/*.pos`, meta-GWAS sumstats, LDREF | `TWAS_results/<tissue>_chr<N>.dat` |
| 2. Run FUSION TWAS | `submit_all_twas_ct.sh` | CT | `TWAS_pos_ct/*.pos`, meta-GWAS sumstats, LDREF | `TWAS_results_ct/TWAS_metaGWAS_<ct>_chr<N>.dat` |
| 2. Single-target prototype | `run_twas.sh` | tissue | (same) | per-chr `.dat` in cwd (whole_blood example) |
| 2. Single-target prototype | `run_twas_ct.sh` | CT | (same) | `TWAS_results_ct/TWAS_metaGWAS_<CT>_chr<N>.dat` |
| 3. QC | `QC_twas_ct.R` | CT | `TWAS_results_ct/*.dat`, `TWAS_pos_ct/*.pos` | `twas_qc_genelist_ct.csv`, `twas_qc_summary_ct.csv` |
| 4. Rebuild path table | `Rebuild_twas_path_table_39.R` | tissue | `TWAS_results/`, `TissueGroups.txt`, heritable-genes filter dirs (N320 / Nall) | `rebuilt/` (gtex_twas + transcripts + `twas_path_table_39_new.csv`) |
| 4. Rebuild path table | `Rebuild_twas_path_table_ct.R` | CT | `TWAS_results_ct/`, `TWAS_pos_ct/`, `twas_qc_genelist_ct.csv` | `rebuilt_ct/` (ct_twas + transcripts + `twas_path_table_ct_new.csv`) |
| 5. Sorted-gene table for FOCUS | `Make_sorted_gene_from_TWAS_results.R` | tissue | `TWAS_results/`, `TissueGroups.txt`, gene-annotation | `FOCUS/sorted_gene_all39Tissue_new.csv` |
| 5. Sorted-gene table for FOCUS | `Make_sorted_gene_from_TWAS_results_ct.R` | CT | `TWAS_results_ct/`, `TWAS_pos_ct/`, `twas_qc_genelist_ct.csv`, gene-annotation | `FOCUS/sorted_gene_all17CT_new.csv` |
| reference | `pos_tissues.txt` | tissue | -- | list of expected tissue pos basenames (39 entries) |

## How to run

All commands assume you are at `${PROJECT_DIR}/TWAS/`
(the scripts hardcode absolute paths under that root).

### Tissue track

```bash
# 1. Build .pos files (TODO: producer script)

# 2. Submit one TWAS job per tissue (each loops chr 1-22).
bash submit_all_twas.sh

# 4. Rebuild the path table for stage-1 PTRS scoring.
Rscript Rebuild_twas_path_table_39.R

# 5. Build the sorted-gene table consumed by FOCUS.
Rscript Make_sorted_gene_from_TWAS_results.R
```

### CT track

```bash
# 1. Build .pos files for the 17 OneK1K cell types.
Rscript Make_Pos_File_celltype.R

# 2. Submit one TWAS job per cell type.
bash submit_all_twas_ct.sh

# 3. QC (HSQ>0, nonzero predictor variance, MODELCV.R2>0).
Rscript QC_twas_ct.R

# 4. Rebuild the CT path table for stage-1 PTRS scoring.
Rscript Rebuild_twas_path_table_ct.R

# 5. Build the sorted-gene table consumed by FOCUS.
Rscript Make_sorted_gene_from_TWAS_results_ct.R
```

The QC step only runs in the CT track. Tissue weights come from the upstream
TCSC heritable-genes filter (`N320` / `Nall` `TranscriptsIn<Tissue>Model_keep.txt`
files), which `Rebuild_twas_path_table_39.R` references directly via
`transcript_keep_file`.

## Inputs expected

- **Meta-analyzed asthma GWAS sumstats**:
  `formatted_meta_analysis.txt` (produced by step 01_meta-analysis)
- **LD reference**: `LDREF/1000G.EUR.` (per-chromosome PLINK bfile prefix)
- **FUSION TWAS toolkit**: `fusion_twas/FUSION.assoc_test.R`
- **GTEx FUSION weights** (tissue track):
  `${REF_DIR}/gtex/weights/v8_320EUR/META_<T>/` and
  `${REF_DIR}/gtex/weights/v8_allEUR_<T>_blup/` (Nall for small tissues)
- **OneK1K celltype weights** (CT track):
  `${REF_DIR}/1k1k_cell_type/celltype_weights_v2/<CT>/`
- **Gene annotation**: `${REF_DIR}/gene_annotation.txt.gz`
  (used by stages 1, 5)
- **TCSC tissue grouping + heritable-genes filter** (tissue track):
  `${REF_DIR}/TCSC/analysis/TissueGroups.txt` and
  `${REF_DIR}/TCSC/weights/heritablegenes/{N320,Nall}/`

## Data inputs (real vs sample candidates)

These are the categories of input data this stage consumes. Sizes are
approximate; "Sample viable?" flags which inputs are small or anonymizable
enough that we could ship a committed example fixture vs which would have to
remain external downloads. The final real-vs-sample choice for each row is
left for the maintainer to fill in.

| Input | Source | Approx. size | Sample viable? | Notes |
|---|---|---|---|---|
| Meta-GWAS sumstats (`formatted_meta_analysis.txt`) | `01_meta-analysis/` output | ~24 MB | **Yes** | Already small; reuse the sample produced by `01_meta-analysis/sample_data/` |
| LD reference (1000G EUR, per-chr PLINK) | external (Hapmap3 / 1000G) | ~500 MB total (~20 MB / chr) | **Maybe** | One chromosome (e.g., chr22) is committable; full set is not |
| FUSION TWAS toolkit | external (gusevlab/fusion_twas) | git clone | N/A | Link to upstream repo only |
| GTEx FUSION weights (39 tissues) | external (GTEx v8 / TCSC release) | ~10 GB | **No** | Document download URL; too large to ship |
| OneK1K celltype weights (17 CTs) | external collaborator | ~5 GB | **No** | Same; document provenance |
| Gene annotation (`gene_annotation.txt.gz`) | external | ~5 MB | **Yes** | Small reference table; committable |
| TCSC tissue grouping (`TissueGroups.txt`) | external | ~KB | **Yes** | Committable |
| TCSC heritable-genes filter (`TranscriptsIn*Model_keep.txt`) | external | ~50 MB total | **Yes** | Per-tissue text files; committable in full or per-tissue |
| Per-CT OneK1K `.fam` files | external | ~MB / CT | **Yes** | Needed only to count `N` for the `.pos` files; an anonymized / shape-only fam file suffices |

## Output schema cheatsheet

- `TWAS_pos{,_ct}/*.pos` — tab-separated, header
  `PANEL  WGT  ID  CHR  P0  P1  N`. One row per gene.
- `TWAS_results{,_ct}/*.dat` — FUSION assoc output (default columns: `PANEL,
  FILE, ID, CHR, P0, P1, HSQ, BEST.GWAS.ID, BEST.GWAS.Z, EQTL.ID, EQTL.R2,
  EQTL.Z, EQTL.GWAS.Z, NSNP, MODEL, MODELCV.R2, MODELCV.PV, TWAS.Z, TWAS.P`).
- `rebuilt/twas_path_table_39_new.csv` — 6 cols: `tissue, gtex_twas,
  transcripts, transcript_keep_file, weights, PANEL`. Quoted.
- `rebuilt_ct/twas_path_table_ct_new.csv` — 5 cols: `tissue, gtex_twas,
  transcripts, weights, PANEL`. Quoted. (No `transcript_keep_file` — QC has
  already been applied.)
- `twas_qc_genelist_ct.csv` — 2 cols: `PANEL, gene_id`. One row per passing
  gene-CT pair.
- `twas_qc_summary_ct.csv` — per-CT counts: `PANEL, total, HSQ_gt0, nzVar,
  R2_gt0, ALL3`.
- `Marginal_alphas_NEW_TWAS_<T>.txt.gz` + `TranscriptsIn<T>Model.txt` — paired
  files, one row per gene, row-aligned (cbind-able). The first carries
  `TWAS.Z`; the second carries the gene ID.
