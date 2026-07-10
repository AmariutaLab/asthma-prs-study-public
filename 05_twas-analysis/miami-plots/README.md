# Miami plots (Figure 3A–B)

Genome-wide TWAS Miami plots of asthma associations for **39 GTEx bulk tissues**
and **17 OneK1K pseudobulked PBMC cell types**. Panel A of the manuscript shows
esophagus mucosa; panel B shows CD4⁺ naïve T cells. The same code also emits a
per-tissue / per-cell-type plot for every panel plus the multi-panel supplement
figures.

Two steps:

1. **`twas_preprocessing.ipynb`** — concatenates the per-(gene, chromosome) FUSION
   TWAS `.dat` outputs into one tidy table per modality, cleans column types, and
   maps tissue / cell-type labels.
2. **`miami_plot.ipynb`** — reads those tables and draws the Miami plots.

## Folder layout

```
miami-plots/
├── twas_preprocessing.ipynb
├── miami_plot.ipynb
├── data/        # inputs + generated tables  (git-ignored)
└── svg/         # output figures             (git-ignored)
```

## Inputs (run on EXPANSE)

The raw inputs are the FUSION TWAS `.dat` files produced upstream in
`05_twas-analysis/` (one per gene-set × chromosome). They are large and live on
EXPANSE; point `data/` at them (symlinks are simplest on the cluster):

```bash
cd 05_twas-analysis/miami-plots
mkdir -p data svg
ln -s /expanse/lustre/projects/ddp412/n5huang/TWAS/TWAS_results    data/TWAS_tissue_results
ln -s /expanse/lustre/projects/ddp412/n5huang/TWAS/TWAS_results_ct data/TWAS_celltype_results
```

`twas_preprocessing.ipynb` globs:
- tissues:    `data/TWAS_tissue_results/twas_*_chr*.dat` (+ `…_chr6.dat.MHC`)
- cell types: `data/TWAS_celltype_results/TWAS_metaGWAS_*_chr*.dat` (+ `…_chr6.dat.MHC`)

so make sure your FUSION output filenames match those patterns.

## Run order

```bash
jupyter nbconvert --to notebook --execute --inplace twas_preprocessing.ipynb
jupyter nbconvert --to notebook --execute --inplace miami_plot.ipynb
```

`twas_preprocessing.ipynb` writes:
- `data/twas_all-tissue-all-chr.tsv`
- `data/twas_all-celltype-all-chr.tsv`

`miami_plot.ipynb` writes to `svg/`:
- `svg/miami_TISSUE_<tissue>.svg`, `svg/miami_CELLTYPE_<celltype>.svg` (one per panel)
- `svg/miami_TISSUE_esophagus_mucosa.svg` (Fig 3A), `svg/miami_CELLTYPE_cd4_naive.svg` (Fig 3B)
- `svg/miami_supplement_TISSUES-PART-<n>.svg`, `svg/miami_CELLTYPES-PART-<n>.svg` (supplement panels)

## Plotting-only (no raw `.dat`)

`miami_plot.ipynb` only needs the two `data/twas_all-*-all-chr.tsv` tables. To plot
on another machine, run step 1 on the server, copy just those two TSVs into `data/`,
and run step 2.

Requires: pandas, numpy, matplotlib. Original notebooks by Michelle Franc Ragsac.
