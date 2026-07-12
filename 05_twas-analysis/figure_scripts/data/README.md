# Input data for the eQTL-weight "cards" figures (Figure 2 area)

These CSVs are the extracted FUSION eQTL weights that the figure builders read
directly. They are bundled here so the figures render from a fresh clone with no
access to the (large, un-versioned) raw `.wgt.RDat` weight archives.

| File | Columns | Focal gene | Source models | Built by | Consumed by |
|---|---|---|---|---|---|
| `gsdmb_blup4.csv` | `tissue, snp, weight` | GSDMB (`ENSG00000073605`, 17q21) | GTEx v8 320EUR, BLUP column | `extract_gtex_weights.R` | `build_eqtl_weight_cards.py` |
| `hla_dpb1_onek1k_blup.csv` | `celltype, snp, weight` | HLA-DPB1 | OneK1K ma-focus, BLUP column | `extract_onek1k_weights.R` | `build_eqtl_weight_cards_1k1k.py` |

## Render the figures (uses these CSVs by default)

```bash
cd 05_twas-analysis/figure_scripts
python build_eqtl_weight_cards.py        # -> eqtl_weight_cards.png
python build_eqtl_weight_cards_1k1k.py   # -> eqtl_weight_cards_1k1k.png
```

Output dir defaults to `$FIG_OUTDIR` (set it to control where the PNGs land).

## Regenerate the CSVs from raw weights

Requires the FUSION `.wgt.RDat` weight files (not in this repo). On the analysis
host they live under `FOCUS/output_39new/v8_320EUR/META_<tissue>/` (GTEx) and
`FOCUS/ma-focus/<CellType>/` (OneK1K). Point the extractors at them:

```bash
# GTEx (GSDMB)
HARTWELL_ROOT=<root> \
WEIGHTS_DIR=<root>/FOCUS/output_39new/v8_320EUR \
OUT_CSV=data/gsdmb_blup4.csv \
Rscript extract_gtex_weights.R

# OneK1K (HLA-DPB1)
HARTWELL_ROOT=<root> \
MAFOCUS_DIR=<root>/FOCUS/ma-focus \
OUT_CSV=data/hla_dpb1_onek1k_blup.csv \
Rscript extract_onek1k_weights.R
```

`METHOD` (default `blup`), `TISSUES` / `CELLTYPES`, `ENSG` / `GENE` are also
overridable via environment variables — see the headers of the two `.R` scripts.
