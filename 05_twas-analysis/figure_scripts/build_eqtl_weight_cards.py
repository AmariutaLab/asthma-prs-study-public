"""
eQTL-weight "cards" (eQTL-weights input panel for the methods / Figure 2 area).

Renders the TWAS eQTL-weight input as a stack of gene-set cards: one card =
one gene, rows = cis-SNPs, columns = tissues, cells = the (scaled) FUSION eQTL
weight. Filled with REAL GTEx weights for GSDMB (17q21 asthma locus); the
trailing cards/columns/rows are shown abstractly (ellipses) to convey
"many genes x many tissues" without clutter.

Visual encoding:
    cell value   eQTL weight scaled to [-1, 1] (per gene); negatives in red
    teal table   light teal grid, darker teal header row, light row-label col
    card stack    two faint cards behind the front card = "more gene sets"

Input (see CONFIG): a long CSV of real weights produced by the companion
    extract_gtex_weights.R  ->  columns: tissue, snp, weight

Outputs (see CONFIG):
    eqtl_weight_cards.png   (300 dpi, transparent)
"""
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

# Editable text in SVG / PDF -- keep glyph data, don't outline to paths
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["pdf.fonttype"] = 42

# =============================================================================
# CONFIG -- edit these paths for your environment
# =============================================================================
ROOT = Path(os.environ.get("HARTWELL_ROOT", "/Users/nancyh/Desktop/hartwell"))
WEIGHTS_CSV = Path(os.environ.get(
    "GSDMB_WEIGHTS_CSV", Path(tempfile.gettempdir()) / "gsdmb_blup4.csv"))
OUTDIR = Path(os.environ.get("FIG_OUTDIR", ROOT / "FOCUS"))
OUTDIR.mkdir(parents=True, exist_ok=True)

GENE = "GSDMB"
TCOLS = ["Lung", "Whole_Blood", "Esophagus_Mucosa", "Thyroid"]   # 4 columns shown
SHORT = {"Lung": "Lung", "Whole_Blood": "Blood",
         "Esophagus_Mucosa": "Esoph.", "Thyroid": "Thyroid"}

if not WEIGHTS_CSV.exists():
    raise SystemExit(
        f"weights CSV not found: {WEIGHTS_CSV}\n"
        f"Generate it first, e.g.:\n"
        f"    Rscript {Path(__file__).with_name('extract_gtex_weights.R')}")

# =============================================================================
# Real weights -> scaled matrix, pick informative SNPs
# =============================================================================
d = pd.read_csv(WEIGHTS_CSV)
M = d.pivot_table(index="snp", columns="tissue", values="weight", aggfunc="first")[TCOLS]
M = M / np.nanmax(np.abs(M.values))                       # scale to [-1, 1]
sel = M.loc[M.var(axis=1).sort_values(ascending=False).head(4).index]   # 4 SNPs
snps = list(sel.index)

col_head = [SHORT[TCOLS[0]], SHORT[TCOLS[1]], "⋯", SHORT[TCOLS[3]]]
row_lab = [snps[0], snps[1], "⋮", snps[3]]


def fmt(v):
    return "NA" if np.isnan(v) else f"{v:+.2f}".replace("+", " ")


def cellval(i, j):
    if row_lab[i] == "⋮":
        return "⋯"
    if col_head[j] == "⋯":
        return "⋮"
    tj = [TCOLS[0], TCOLS[1], None, TCOLS[3]][j]
    return fmt(sel.loc[row_lab[i], tj])


# =============================================================================
# Draw the card stack + table
# =============================================================================
TEAL_HDR, TEAL_CELL, TEAL_LAB = "#cfe8e4", "#eaf6f4", "#dcefeb"
EDGE, TXT, HDRTXT = "#9ec9c3", "#37474f", "#2f6f68"

fig = plt.figure(figsize=(3.4, 2.5), dpi=300)
fig.patch.set_alpha(0)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

for off in [(10, -8), (5, -4)]:                          # faint back cards
    ax.add_patch(FancyBboxPatch((14 + off[0], 14 + off[1]), 80, 74,
                 boxstyle="round,pad=0,rounding_size=3", linewidth=1.1,
                 edgecolor=EDGE, facecolor="#f1faf8", zorder=1))
ax.add_patch(FancyBboxPatch((14, 14), 80, 74, boxstyle="round,pad=0,rounding_size=3",
             linewidth=1.3, edgecolor="#6fb3aa", facecolor="white", zorder=5))
ax.text(18, 92, "✕", fontsize=8, color="#8aa", zorder=6, va="center")
ax.text(22, 92, GENE, fontsize=10, fontweight="bold", fontstyle="italic",
        color=HDRTXT, zorder=6, va="center")
ax.text(7.5, 52, "Gene set A", rotation=90, fontsize=7.5, color="#6b8e89",
        ha="center", va="center", zorder=6)

x0, y0, x1, y1 = 19, 19, 90, 84
nC, nR = 5, 5
cw = (x1 - x0) / nC
rh = (y1 - y0) / nR
cx = lambda c: x0 + c * cw
cy = lambda r: y1 - (r + 1) * rh

for r in range(nR):
    for c in range(nC):
        fc = ("none" if (r == 0 and c == 0) else TEAL_HDR if r == 0
              else TEAL_LAB if c == 0 else TEAL_CELL)
        ax.add_patch(Rectangle((cx(c), cy(r)), cw, rh, facecolor=fc,
                     edgecolor="white", linewidth=1.1, zorder=6))
        if r == 0 and c == 0:
            continue
        if r == 0:                                       # column headers
            ax.text(cx(c) + cw / 2, cy(r) + rh / 2, col_head[c - 1], ha="center",
                    va="center", fontsize=6.6, fontweight="bold", color=HDRTXT, zorder=7)
        elif c == 0:                                     # row labels
            t = row_lab[r - 1]
            ax.text(cx(c) + cw / 2, cy(r) + rh / 2, t, ha="center", va="center",
                    fontsize=5.4 if t.startswith("rs") else 7, fontweight="bold",
                    color="#5a7d78", zorder=7)
        else:                                            # values
            v = cellval(r - 1, c - 1)
            col = ("#b94b4b" if (v not in ("NA", "⋯", "⋮") and v.strip().startswith("-"))
                   else "#9aa" if v == "NA" else TXT)
            ax.text(cx(c) + cw / 2, cy(r) + rh / 2, v, ha="center", va="center",
                    fontsize=6.5, color=col, zorder=7)

png = OUTDIR / "eqtl_weight_cards.png"
fig.savefig(png, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.02)
print(f"saved {png}")
print(sel.round(2).to_string())
