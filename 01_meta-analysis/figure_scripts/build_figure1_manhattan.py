"""
Figure 1, Panel A -- GBMI+TAGC asthma meta-analysis Manhattan plot.

Renders the genome-wide Manhattan for the fixed-effects IVW meta-analysis:
alternating silver/grey autosomes, a dashed 5e-8 significance line, index-variant
loci highlighted in orange-red, and the top-3 nearest-gene symbols labelled at
each independent locus (rsIDs stripped from the labels).

Pipeline:
    1. Load the formatted meta-analysis sumstats; compute -log10(p).
    2. Merge dbSNP (b151 / GRCh37) annotations to recover GENEINFO gene symbols.
    3. Compute a running genomic "Position" with padding between chromosomes.
    4. Flag index variants from the PLINK clumping result (r2 0.1 / 1 Mb).
    5. Plot all significant SNPs + a 10% sample of the rest (seed 9550);
       highlight index-variant loci and annotate their top-3 genes.

Inputs (all aggregate GWAS / annotation data -- NOT individual-level; see CONFIG):
    <FIG1_DATA_ROOT>/summary-stats/metaanalysis_merged-output.formatted.pval.txt
    <FIG1_DATA_ROOT>/human_9606_b151_GRCh37p13.filt-metanaalysis-var.txt   (dbSNP)
    <FIG1_DATA_ROOT>/plink-clumping-results/
        VERBOSE_metaanalysis-index-var.allchr.clump_r2-0.1_1000-kb.clumped

Outputs (see CONFIG -> FIG_OUTDIR, default = ../figures):
    figure1_manhattan_labeled.svg   (editable text, dpi 1200)
    figure1_manhattan_labeled.png   (300 dpi preview)
    figure1_manhattan_annotated-INDEXVAR-ONLY.tsv   (index-variant table; reused by Fig 2)

NB: this is the *Figure 1* Manhattan. The Figure 2 methods-overview Manhattan is a
separate script, `build_manhattan_plot.py`.
"""
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

# Keep text as text (not outlines) in SVG / PDF exports
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"


# =============================================================================
# CONFIG -- edit paths for your environment (env vars override the defaults)
# =============================================================================
DATA_ROOT = Path(os.environ.get(
    "FIG1_DATA_ROOT",
    "/Users/nancyh/Desktop/hartwell-asthma-main/local-computer-files/figure-1/data"))

META_SUMSTATS = Path(os.environ.get(
    "META_SUMSTATS", DATA_ROOT / "summary-stats" /
    "metaanalysis_merged-output.formatted.pval.txt"))

DBSNP_ANNOT = Path(os.environ.get(
    "DBSNP_ANNOT", DATA_ROOT / "human_9606_b151_GRCh37p13.filt-metanaalysis-var.txt"))

CLUMP_FILE = Path(os.environ.get(
    "CLUMP_FILE", DATA_ROOT / "plink-clumping-results" /
    "VERBOSE_metaanalysis-index-var.allchr.clump_r2-0.1_1000-kb.clumped"))

FIG_OUTDIR = Path(os.environ.get(
    "FIG_OUTDIR", Path(__file__).resolve().parents[1] / "figures"))
FIG_OUTDIR.mkdir(parents=True, exist_ok=True)

SIGNIFICANCE = float(os.environ.get("SIGNIFICANCE", "5e-8"))
CHROM_PADDING = 4e7    # blank space between chromosomes on the x-axis
IVAR_FLANK = 5e5       # +/- window (bp) grouping points around an index variant
SAMPLE_SEED = 9550     # non-significant point subsample seed


# =============================================================================
# Step 1 -- load + annotate meta-analysis sumstats
# =============================================================================
for f in (META_SUMSTATS, DBSNP_ANNOT, CLUMP_FILE):
    if not f.exists():
        raise FileNotFoundError(f"Required input not found: {f}")

df_meta = pd.read_csv(META_SUMSTATS, delimiter="\t")
df_meta["-log10(p-value)"] = df_meta["pval"].apply(lambda r: -1 * np.log10(float(r)))

df_dbsnp = pd.read_csv(DBSNP_ANNOT, delimiter="\t")
df = pd.merge(df_meta, df_dbsnp, left_on="SNP", right_on="ID").sort_values("pval")

colnames = ["SNP", "#CHROM", "POS", "A2", "A1", "REF", "ALT",
            "pval", "-log10(p-value)", "BETA", "SE", "Z", "INFO"]
df = df[colnames].rename(columns={"pval": "p-value"})
df = df.loc[df["#CHROM"] != "X"]
df["#CHROM"] = pd.to_numeric(df["#CHROM"])


# =============================================================================
# Step 2 -- running genomic position (with inter-chromosome padding)
# =============================================================================
chrom_lengths = df.sort_values(["#CHROM", "POS"]).groupby("#CHROM")["POS"].max() + CHROM_PADDING
chrom_lengths[22] = chrom_lengths[22] - CHROM_PADDING   # no padding after the last chromosome
chrom_offsets = chrom_lengths.cumsum().shift(1, fill_value=0)
df["Position"] = df["POS"] + df["#CHROM"].map(chrom_offsets)
chrom_midpoints = df.groupby("#CHROM")["Position"].mean()


# =============================================================================
# Step 3 -- gene symbols from the dbSNP GENEINFO field
# =============================================================================
df["GENEINFO"] = df.apply(
    lambda r: re.search(r"GENEINFO=(.*?);", r["INFO"]).group(1)
    if "GENEINFO" in r["INFO"] else r["SNP"], axis=1)
df["GENE_SYMBOLS"] = df["GENEINFO"].apply(lambda r: [g.split(":")[0] for g in r.split("|")])


# =============================================================================
# Step 4 -- flag index variants from the clumping result
# =============================================================================
header_col = ["CHR", "F", "SNP", "BP", "P", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001"]
ignore_col = {"CHR", "F", "BP", "P"}
ivar_rsids = set()
with open(CLUMP_FILE) as fh:
    take = False
    for line in fh:
        if "CHR" in line:
            take = True
            continue
        if take and line.strip():
            rec = {c: v for c, v in zip(header_col, line.strip().split()) if c not in ignore_col}
            ivar_rsids.add(rec["SNP"])
            take = False
df["INDEX_VAR"] = df["SNP"].isin(ivar_rsids)
df = df.sort_values(["#CHROM", "POS"]).reset_index(drop=True)
print(f"Meta-analysis: {df.shape[0]:,} annotated SNPs; {df['INDEX_VAR'].sum()} index variants")


# =============================================================================
# Step 5 -- plotting frame: all significant + 10% of the rest
# =============================================================================
sig = df[df["p-value"] <= SIGNIFICANCE]
nonsig = df[df["p-value"] > SIGNIFICANCE].sample(frac=0.1, random_state=SAMPLE_SEED)
df_plot = pd.concat([sig, nonsig]).sort_values("Position")
print(f"Plotting {df_plot.shape[0]:,} points ({sig.shape[0]:,} significant)")

chrom_colors = ["silver", "dimgray"]
cmap = {c: chrom_colors[i % 2] for i, c in enumerate(df_plot["#CHROM"].unique())}
df_plot["CHROM_COLOR"] = df_plot["#CHROM"].map(cmap)

# Group points around each locus's index variant (top 10 per chromosome) and
# collect the top-3 non-rsID gene symbols for the locus label.
df_ivar = (df_plot.loc[df_plot["INDEX_VAR"]]
           .sort_values(["#CHROM", "p-value"]).groupby("#CHROM").head(10))
labels, seen = [], set()
for _, row in df_ivar.iterrows():
    pts = df_plot.loc[(df_plot["Position"] > row["Position"] - IVAR_FLANK) &
                      (df_plot["Position"] < row["Position"] + IVAR_FLANK)]
    if any(p in seen for p in pts.index):
        continue
    seen.update(pts.index)
    genes, seen_g = [], set()
    for g in pts.loc[pts["p-value"] < SIGNIFICANCE].sort_values("p-value")["GENE_SYMBOLS"].sum():
        if g not in seen_g:
            seen_g.add(g)
            if "rs" not in g:
                genes.append(g)
    labels.append({"Position": pts["Position"].mean(),
                   "-log10(p-value)": pts["-log10(p-value)"].max(),
                   "GENE_SYMBOLS": genes[:3]})
ivar_points = df_plot.loc[list(seen)]
ivar_labels = pd.DataFrame(labels)


# =============================================================================
# Step 6 -- render
# =============================================================================
var_point_s, xaxis_padding, ivar_color = 4, 5e7, "orangered"
fig, ax = plt.subplots(1, 1, figsize=(10, 4))
ax.scatter(df_plot["Position"], df_plot["-log10(p-value)"],
           c=df_plot["CHROM_COLOR"], s=var_point_s, alpha=.6, rasterized=True)
ax.scatter(ivar_points["Position"], ivar_points["-log10(p-value)"],
           c=ivar_color, s=var_point_s + 1, alpha=.7, rasterized=True)
ax.axhline(-np.log10(SIGNIFICANCE), color="red", linestyle="--", label="p-value = 5e-8")

for num, row in ivar_labels.iterrows():
    ax.text(row["Position"] - 2e7, row["-log10(p-value)"] + 10,
            ", ".join(row["GENE_SYMBOLS"]), fontstyle="italic", fontsize="x-small",
            rotation=90, gid=f"{num}")   # gid groups each locus for SVG editing

ax.set_xticks(chrom_midpoints)
ax.set_xticklabels(chrom_midpoints.index.astype(int))
ax.set_xlim(df_plot["Position"].min() - xaxis_padding, df_plot["Position"].max() + xaxis_padding)
ax.set_xlabel("Chromosome")
ax.set_ylim(df_plot["-log10(p-value)"].min() - 2, df_plot["-log10(p-value)"].max() + 10)
ax.set_ylabel(r"$-\log_{10}(p-value)$")
plt.tight_layout()

out_svg = FIG_OUTDIR / "figure1_manhattan_labeled.svg"
out_png = FIG_OUTDIR / "figure1_manhattan_labeled.png"
plt.savefig(out_svg, format="svg", dpi=1200)
plt.savefig(out_png, dpi=300)

# Index-variant table (aggregate; reused when assembling the Figure 2 overview)
ivar_tbl = df.loc[df["INDEX_VAR"]].copy()
ivar_tbl["LOCUSZOOM_COORDINATE"] = ivar_tbl.apply(lambda r: f"chr{r['#CHROM']}:{r['POS']}", axis=1)
out_tsv = FIG_OUTDIR / "figure1_manhattan_annotated-INDEXVAR-ONLY.tsv"
ivar_tbl.to_csv(out_tsv, sep="\t", index=False)

print(f"Saved: {out_svg}\n       {out_png}\n       {out_tsv}")
