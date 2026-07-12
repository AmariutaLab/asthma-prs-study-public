"""
Single-cell eQTL-weight "cards" (OneK1K eQTL-weights input panel, methods /
Figure 2 area) -- the cell-type counterpart of build_eqtl_weight_cards.py.

Renders the TWAS eQTL-weight input as a stack of gene-set cards: one card =
one gene, rows = cis-SNPs, columns = OneK1K cell types, cells = the (scaled)
FUSION eQTL weight. Filled with REAL OneK1K weights for HLA-DPB1 (the cell-type
panel focal gene); trailing cards/columns/rows are shown abstractly (ellipses)
to convey "many genes x many cell types" without clutter.

Visual encoding:
    cell value   eQTL weight scaled to [-1, 1] (per gene); negatives in red
    teal table   light teal grid, darker teal header row, light row-label col
    card stack    two faint cards behind the front card = "more gene sets"

Input (see CONFIG): a long CSV of real weights produced by the companion
    extract_onek1k_weights.R  ->  columns: celltype, snp, weight

Outputs (see CONFIG):
    eqtl_weight_cards_1k1k.png   (300 dpi, transparent)
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
# Default to the extracted CSV bundled next to this script (data/), so a fresh
# clone renders out of the box; override with ONEK1K_WEIGHTS_CSV to point at a
# freshly re-extracted file (see extract_onek1k_weights.R / data/README.md).
WEIGHTS_CSV = Path(os.environ.get(
    "ONEK1K_WEIGHTS_CSV", Path(__file__).with_name("data") / "hla_dpb1_onek1k_blup.csv"))
OUTDIR = Path(os.environ.get("FIG_OUTDIR", ROOT / "FOCUS"))
OUTDIR.mkdir(parents=True, exist_ok=True)

GENE = "HLA-DPB1"
# 4 OneK1K cell-type columns shown (CamelCase dir names == values in the CSV)
CCOLS = ["B_memory", "CD4_CTL", "CD4_TCM", "NK"]
SHORT = {"B_memory": "B mem", "CD4_CTL": "CD4 CTL",
         "CD4_TCM": "CD4 TCM", "NK": "NK"}

if not WEIGHTS_CSV.exists():
    raise SystemExit(
        f"weights CSV not found: {WEIGHTS_CSV}\n"
        f"Generate it first, e.g.:\n"
        f"    Rscript {Path(__file__).with_name('extract_onek1k_weights.R')}")

# =============================================================================
# Real weights -> scaled matrix, pick informative SNPs
# =============================================================================
d = pd.read_csv(WEIGHTS_CSV)
M = d.pivot_table(index="snp", columns="celltype", values="weight", aggfunc="first")[CCOLS]
M = M / np.nanmax(np.abs(M.values))                       # scale to [-1, 1]
# HLA-DPB1 has a real eQTL model only in the strongest cell type (the others
# are ~0), so the top-variance SNPs are LD-tied at that column's max and give
# near-identical rows. Instead pick 4 SNPs that SPAN the signed weight range of
# the lead (max-variance) cell type, so the displayed rows are distinct and
# representative (a strong +, a mid, and a negative).
lead = M.var(axis=0).idxmax()
order = M[lead].sort_values(ascending=False)
pos = [0, len(order) // 3, 2 * len(order) // 3, len(order) - 1]
sel = M.loc[order.index[pos]]
snps = list(sel.index)

col_head = [SHORT[CCOLS[0]], SHORT[CCOLS[1]], "⋯", SHORT[CCOLS[3]]]
row_lab = [snps[0], snps[1], "⋮", snps[3]]


def fmt(v):
    return "NA" if np.isnan(v) else f"{v:+.2f}".replace("+", " ")


def cellval(i, j):
    if row_lab[i] == "⋮":
        return "⋯"
    if col_head[j] == "⋯":
        return "⋮"
    cj = [CCOLS[0], CCOLS[1], None, CCOLS[3]][j]
    return fmt(sel.loc[row_lab[i], cj])


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
# in-card "✕ HLA-DPB1" title + "Gene set A" side strip intentionally omitted so
# annotations can be added by hand downstream.

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

png = OUTDIR / "eqtl_weight_cards_1k1k.png"
fig.savefig(png, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.02)
print(f"saved {png}")
print(sel.round(2).to_string())
