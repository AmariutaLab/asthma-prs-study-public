# Meta-Analysis Pipeline

Fixed-effects inverse-variance weighted (IVW) meta-analysis combining [TAGC](https://www.ebi.ac.uk/gwas/publications/29273806) and [GBMI](https://www.globalbiobankmeta.org/) asthma summary statistics, with downstream formatting for TWAS and PRS analyses.

## Pipeline Overview

```
TAGC sumstats (.tsv)  ──┐                                                   
                        ├─ Stage 1 ─► Standardized (A1, A2, SNP, N, Z)      
GBMI sumstats (.txt.gz) ┘                                                   
        │                                                                    
        ├─ Stage 2 ─► IVW meta-analysis (common SNPs only)                  
        │              output: meta_analysis_results_common.csv              
        │                                                                    
        ├─ Stage 3 ─► Outer-join to recover study-specific SNPs             
        │              output: meta_analysis_output.csv                      
        │                                                                    
        └─ Stage 4 ─► Format for downstream (TWAS, PRS)                     
                       output: formatted_meta_analysis.txt                   
                               meta_analysis_combineStat.txt                 
                               meta_analysis_p_complete.txt.gz               
```

## Allele Convention

Consistent across all stages:

| Column | Meaning |
|--------|---------|
| **A1** | ALT / alternate allele (effect allele) |
| **A2** | REF / reference allele |

Beta values are always relative to A1.

## Scripts

| Script | Stage | Description |
|--------|-------|-------------|
| `01_format_GBMI.R` | 1 | Standardize GBMI columns to `A1, A2, SNP, N, Z` |
| `01_format_TAGC.R` | 1 | Standardize TAGC columns to `A1, A2, SNP, N, Z` |
| `02_meta_common.R` | 2 | IVW fixed-effects meta-analysis on common SNPs (uses `metafor::rma`) |
| `03_merge_all_snps.py` | 3 | Outer-join to include GBMI-only and TAGC-only SNPs |
| `04_format_downstream.py` | 4 | Assign alleles, compute Z/p-values, format for TWAS and PRS |
| `run_test.py` | Test | Sample 100 GBMI + 50 TAGC SNPs, run pipeline, compare with full results |

## Output Files

| File | Columns | Rows | Used By |
|------|---------|------|---------|
| `meta_analysis_results_common.csv` | chr, rsid, beta/se (per study + combined) | ~1.6M (common SNPs) | Stage 3 |
| `meta_analysis_output.csv` | chr, rsid, combined_beta/se, gbmi_beta/se, tagc_beta/se | ~28M (all SNPs) | Stage 4 |
| `meta_analysis_combineStat.txt` | SNP, BETA, SE | ~28M | (archival intermediate; not consumed by downstream scripts) |
| `formatted_meta_analysis.txt` | SNP, Z, A2, A1 | ~24M | TWAS / FOCUS |
| `meta_analysis_p_complete.txt.gz` | SNP, A1, A2, Z, #CHR, POS, BETA, SE, p_value, -log10(P-Value) | ~24M | Clumping / PRS |

## Figures

Standalone plotting scripts live in [`figure_scripts/`](figure_scripts/). Each
reads its paths from environment variables (defaults point at the original lab
layout) so it runs unchanged on a laptop or the server.

| Script | Output | Description |
|--------|--------|-------------|
| `build_figure1_manhattan.py` | `figure1_manhattan_labeled.svg` (raw panel; regenerated locally, not committed) | **Figure 1, Panel A** — GBMI+TAGC meta-analysis Manhattan: alternating-grey autosomes, dashed 5e-8 line, index-variant loci highlighted in orange-red with top nearest-gene labels. Reads the formatted sumstats + dbSNP (b151) annotations + the PLINK clumping result. The script emits the panel with **vertical** gene labels; the **published** Panel A was then hand-curated (labels rotated horizontal, grouped, leader lines to peaks) in Affinity Designer. The assembled Figure 1 (Panel A + Panel B) is committed as [`figures/figure1_combined.png`](figures/figure1_combined.png) / `.pdf`, with the editable Affinity source `figures/figure1_combined.af`. |
| `build_manhattan_plot.py` | [`figures/figure2_manhattan_labeled.png`](figures/figure2_manhattan_labeled.png) / `.svg` | GBMI+TAGC EUR meta-analysis Manhattan (Figure 2 overview): genome-wide-significant peaks + curated asthma-locus gene labels over 23,136,851 SNPs (23,145 genome-wide significant). Reads the full sumstats (`meta_analysis_p_complete.txt.gz`) + nearest-gene labels; committed under [`figures/`](figures/). |
| `build_ld_reference_heatmap.py` | [`figures/ld_reference_1kg.png`](figures/ld_reference_1kg.png) / `.pdf` | Lower-triangle r² LD heatmap over the top-N cis-SNPs (ranked by max \|FUSION eQTL weight\| across GTEx tissues) at the *GSDMB* / 17q21 locus, using the 1000G EUR reference panel (`1000G.EUR.17`). Committed under [`figures/`](figures/). |

### `build_figure1_manhattan.py` — environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `FIG1_DATA_ROOT` | `…/hartwell-asthma-main/local-computer-files/figure-1/data` | Base dir the input defaults hang off |
| `META_SUMSTATS` | `$FIG1_DATA_ROOT/summary-stats/metaanalysis_merged-output.formatted.pval.txt` | Formatted meta-analysis sumstats (aggregate GWAS) |
| `DBSNP_ANNOT` | `$FIG1_DATA_ROOT/human_9606_b151_GRCh37p13.filt-metanaalysis-var.txt` | dbSNP b151/GRCh37 annotations (for GENEINFO gene symbols) |
| `CLUMP_FILE` | `$FIG1_DATA_ROOT/plink-clumping-results/VERBOSE_metaanalysis-index-var.allchr.clump_r2-0.1_1000-kb.clumped` | PLINK clumping result (r² 0.1 / 1 Mb) → index variants |
| `FIG_OUTDIR` | `../figures` | Where the SVG/PNG (+ index-variant TSV) are written |

All inputs are aggregate GWAS/annotation data (not individual-level); they are
large and live on the server, so they are not committed. Requires: `pandas`,
`numpy`, `matplotlib`.

### `build_ld_reference_heatmap.py` — environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `HARTWELL_ROOT` | `/Users/nancyh/Desktop/hartwell` | Base dir the other defaults hang off |
| `GSDMB_WEIGHTS_CSV` | `$TMPDIR/gsdmb_blup4.csv` | Long CSV (`tissue, snp, weight`) from `05_twas-analysis/figure_scripts/extract_gtex_weights.R` |
| `LDREF_BFILE` | `$HARTWELL_ROOT/TWAS/LDREF/1000G.EUR.17` | FUSION 1000G EUR LD reference prefix (chr17) |
| `PLINK_BIN` | `$HARTWELL_ROOT/plink` | PLINK 1.9 binary |
| `FIG_OUTDIR` | `$HARTWELL_ROOT/FOCUS` | Where the PNG/PDF are written |
| `LD_TOP_N` | `26` | Number of top-\|weight\| SNPs in the heatmap |
| `LOCUS_GENE` / `LOCUS_REGION` | `GSDMB` / `17q21` | Labels in the plot title |

Pipeline: rank cis-SNPs by max \|weight\| → `plink --r square` + `--write-snplist`
→ mask the upper triangle → matplotlib r² heatmap (mean off-diagonal r² ≈ 0.27
for the committed GSDMB panel). To regenerate the committed figure:

```bash
GSDMB_WEIGHTS_CSV=path/to/gsdmb_blup4.csv \
LDREF_BFILE=path/to/1000G.EUR.17 PLINK_BIN=$(which plink) \
FIG_OUTDIR=01_meta-analysis/figures \
python 01_meta-analysis/figure_scripts/build_ld_reference_heatmap.py
```

Requires: `pandas`, `numpy`, `matplotlib`, plus PLINK 1.9 on `PLINK_BIN`.

## Testing

```bash
cd 01_meta-analysis
python run_test.py
```

`run_test.py` has two modes, auto-selected based on what's available locally:

- **REAL mode** — when the full GBMI + TAGC sumstats are present on disk, it samples 100 GBMI + 50 TAGC SNPs (30 guaranteed common, seed=42), runs Stages 2–4 on the sample, and compares every output against the corresponding rows in the real full results.
- **SAMPLE mode** — when the full sumstats are absent, it loads the committed sample inputs from `sample_data/`, runs Stages 2–4 on them, and compares against the committed expected outputs in `sample_data/expected_outputs/`. This lets anyone clone the repo and run the test without the ~1 GB source files.

All numeric columns must match to within 1e-6 in either mode.

### `sample_data/` layout

```
sample_data/
├── generate_sample_data.py          # one-shot generator (regenerates from full sumstats)
├── sample_gbmi.txt.gz               # 100 rows from GBMI EUR (30 common, 70 GBMI-only)
├── sample_tagc.tsv                  # 50 rows from TAGC EUR (30 common, 20 TAGC-only)
└── expected_outputs/
    ├── meta_analysis_results_common.csv
    ├── meta_analysis_output.csv
    ├── meta_analysis_combineStat.txt
    ├── formatted_meta_analysis.txt
    └── meta_analysis_p_complete.txt.gz
```

To regenerate `sample_data/` from a fresh draw of the full sumstats:

```bash
cd 01_meta-analysis/sample_data
DATA_ROOT=/path/to/parent/of/GBMI_and_TAGC python generate_sample_data.py
```

## Environment setup

Two runtimes: **R** (Stages 1–2) and **Python** (Stages 3–4 + `run_test.py`).

```r
# R (the study used R 4.3.1)
install.packages(c("data.table", "dplyr", "metafor", "tidyr", "purrr", "progress"))
# metafor 4.8-0 was used
```

```bash
# Python
conda create -n meta-analysis -c conda-forge python=3.11 pandas numpy scipy
conda activate meta-analysis
```

**Dependencies at a glance** — R: `data.table`, `dplyr`, `metafor`, `tidyr`,
`purrr`, `progress`; Python: `pandas`, `numpy`, `scipy`.
