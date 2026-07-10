"""
European-ancestry meta-analysis Manhattan plot (Meta-GWAS panel for the
methods / Figure 2 overview area).

Reproduces the published GBMI+TAGC asthma meta-analysis Manhattan: orange
genome-wide-significant peaks on an alternating-grey autosomal background,
a dashed 5e-8 significance line, and ONLY the curated asthma-locus gene
labels, each anchored at its real locus (max -log10 p among the group's
nearest-gene-annotated SNPs; chromosome-peak fallback for genes absent from
the annotation, e.g. RAD51B).

Inputs (see CONFIG block):
    Clumping/meta_analysis_p_complete.txt.gz   full sumstats with #CHR, POS, p_value
    meta_analysis/significant_snps_with_genes.csv  nearest-gene labels for sig SNPs
      NOTE: in that CSV the columns are mislabeled -- `NearestGeneID` holds the
      readable HGNC symbol (e.g. GSDMB) and `NearestGeneSymbol` holds the ENSG id.

Outputs (see CONFIG block):
    figure2_manhattan_labeled.png  (300 dpi)
    figure2_manhattan_labeled.svg  (vector, editable text)

The first run loads ~23M SNPs (~1 GB gz) and caches the positioned plotting
points to a temp pickle so subsequent styling runs are fast.
"""
import os
import tempfile
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Editable text in SVG / PDF -- keep glyph data, don't outline to paths
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# =============================================================================
# CONFIG -- edit these paths for your environment (local laptop vs server)
# =============================================================================
ROOT = Path(os.environ.get("HARTWELL_ROOT", "/Users/nancyh/Desktop/hartwell"))
P_COMPLETE = Path(os.environ.get(
    "META_P_COMPLETE", ROOT / "Clumping" / "meta_analysis_p_complete.txt.gz"))
SIG_GENES = Path(os.environ.get(
    "SIG_SNPS_GENES", ROOT / "meta_analysis" / "significant_snps_with_genes.csv"))
OUTDIR = Path(os.environ.get("FIG_OUTDIR", ROOT / "FOCUS"))
OUTDIR.mkdir(parents=True, exist_ok=True)
CACHE = Path(os.environ.get(
    "MANHATTAN_CACHE", Path(tempfile.gettempdir()) / "manhattan_cache.pkl"))

significance = 5e-8
SIG_Y = -np.log10(significance)
CHROM_PADDING = 4e7          # gap between chromosomes on the running x-axis
NONSIG_FRAC = 0.06           # fraction of non-significant points kept (speed)
RANDOM_STATE = 9550

# =============================================================================
# Load + position (cached)
# =============================================================================
if CACHE.exists():
    with open(CACHE, "rb") as fh:
        C = pickle.load(fh)
    print(f"loaded cache: {CACHE}")
else:
    df = pd.read_csv(P_COMPLETE, sep="\t").rename(
        columns={"#CHR": "#CHROM", "p_value": "p-value",
                 "-log10(P-Value)": "-log10(p-value)"})
    df["#CHROM"] = pd.to_numeric(df["#CHROM"], errors="coerce")
    df = df.dropna(subset=["#CHROM", "POS", "p-value"])
    df["#CHROM"] = df["#CHROM"].astype(int)
    df = df.loc[df["#CHROM"].between(1, 22)].copy()      # autosomes only
    df["-log10(p-value)"] = -np.log10(df["p-value"].astype(float))

    clen = df.sort_values(["#CHROM", "POS"]).groupby("#CHROM")["POS"].max() + CHROM_PADDING
    clen[clen.index.max()] -= CHROM_PADDING              # no pad after last chrom
    offsets = clen.cumsum().shift(1, fill_value=0)
    df["Position"] = df["POS"] + df["#CHROM"].map(offsets)
    mids = df.groupby("#CHROM")["Position"].mean()

    sig = df[df["p-value"] <= significance]
    nons = df[df["p-value"] > significance].sample(frac=NONSIG_FRAC, random_state=RANDOM_STATE)
    keep = ["SNP", "#CHROM", "POS", "Position", "-log10(p-value)", "p-value"]
    C = {"plot": pd.concat([sig, nons])[keep].sort_values("Position"),
         "offsets": offsets.to_dict(), "mids": mids,
         "sig_n": len(sig), "tot": len(df)}
    with open(CACHE, "wb") as fh:
        pickle.dump(C, fh)
    print(f"built cache: {C['tot']:,} SNPs, {C['sig_n']:,} significant")

dfp = C["plot"]
offsets = C["offsets"]
mids = C["mids"]

# =============================================================================
# Gene anchors from the nearest-gene annotation
# =============================================================================
sg = pd.read_csv(SIG_GENES)
sg["#CHROM"] = pd.to_numeric(sg["#CHR"], errors="coerce")
sg = sg.dropna(subset=["#CHROM", "POS"])
sg["#CHROM"] = sg["#CHROM"].astype(int)
sg["nl10"] = pd.to_numeric(sg["-log10(P-Value)"], errors="coerce")
sg["Position"] = sg["POS"] + sg["#CHROM"].map(offsets)
sg["SYM"] = sg["NearestGeneID"].astype(str).str.upper()   # see docstring note

CHR_FALLBACK = {"RAD51B": 14}    # genes absent from nearest-gene annot -> chrom peak


def anchor(genes):
    """Best (Position, peak height) among the given gene symbols; fall back to
    the strongest peak on a known chromosome for genes missing from the annot."""
    sub = sg[sg["SYM"].isin([g.upper() for g in genes])]
    if not sub.empty:
        r = sub.loc[sub["nl10"].idxmax()]
        return float(r["Position"]), float(r["nl10"])
    for g in genes:
        ch = CHR_FALLBACK.get(g.upper())
        if ch:
            cc = dfp[dfp["#CHROM"] == ch]
            r = cc.loc[cc["-log10(p-value)"].idxmax()]
            return float(r["Position"]), float(r["-log10(p-value)"])
    return None


# Curated label groups from the published figure: (genes, label_y, x_nudge_Mb)
GROUPS = [
    (["PEX14"], 80, 0),
    (["TNFSF4", "IL6R", "ADORA1"], 60, 0),
    (["IL1RL1", "IL18R1"], 104, 0),
    (["D2HGDH"], 58, 0),
    (["LPP"], 40, 0),
    (["TLR1"], 50, 0),
    (["KIAA1109", "IL2"], 86, 0),
    (["TSLP", "WDR36", "IL13", "NDFIP1"], 124, -3),
    (["HLA-DQA1", "HLA-DQB1", "HLA-DQB1-AS1", "MICA"], 146, 0),
    (["BACH2"], 56, 6),
    (["IL6", "CDHR3"], 36, 0),
    (["IL33", "TPD52L3"], 122, 0),
    (["EMSY"], 56, 0),
    (["STAT6"], 46, 0),
    (["RAD51B"], 40, 0),
    (["SMAD3", "RORA"], 70, -4),
    (["CLEC16A", "IL4R", "IL21R"], 102, 0),
    (["GSDMB", "ZPBP2", "STAT5B", "STAT5A"], 134, 3),
    (["IL2RB"], 40, 0),
]

# =============================================================================
# Plot
# =============================================================================
ORANGE = "#E8552D"
GREY_LO, GREY_HI = "#d6d6d6", "#a3a3a3"
chroms = sorted(dfp["#CHROM"].unique())
gmap = {c: (GREY_LO if i % 2 == 0 else GREY_HI) for i, c in enumerate(chroms)}
issig = dfp["p-value"].values <= significance
colors = np.where(issig, ORANGE, dfp["#CHROM"].map(gmap).values)

fig, ax = plt.subplots(figsize=(10, 5.2))
ax.scatter(dfp["Position"].values, dfp["-log10(p-value)"].values,
           c=colors, s=4, alpha=.75, rasterized=True, linewidths=0)
ax.axhline(SIG_Y, color="#7a1f1f", linestyle="--", linewidth=1.1, dashes=(6, 4))

for genes, ly, nx in GROUPS:
    a = anchor(genes)
    if a is None:
        print("NOT FOUND:", genes)
        continue
    px, ph = a
    lx = px + nx * 1e6
    ax.plot([px, lx], [ph + 2, ly - 2.5], color="#9aa0a6", lw=0.6, zorder=2)
    ax.text(lx, ly, "\n".join(genes), fontsize=7.0, fontstyle="italic",
            ha="center", va="bottom", color="#222", linespacing=1.05, zorder=6)

ax.set_xticks([mids[c] for c in chroms])
ax.set_xticklabels([int(c) for c in chroms], fontsize=9)
ax.set_xlim(dfp["Position"].min() - 3e7, dfp["Position"].max() + 3e7)
ax.set_ylim(-2, 152)
ax.set_xlabel("Chromosome", fontsize=12)
ax.set_ylabel(r"$-\log_{10}(p-value)$", fontsize=12)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.tick_params(length=3)
plt.tight_layout(rect=[0, 0.06, 1, 1])
fig.text(0.5, 0.012, "European-ancestry Meta-analysis", ha="center",
         fontweight="bold", fontsize=14)

png = OUTDIR / "figure2_manhattan_labeled.png"
svg = OUTDIR / "figure2_manhattan_labeled.svg"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(svg, bbox_inches="tight")
print(f"saved {png}\n      {svg}")
