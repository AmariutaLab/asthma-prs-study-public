"""
Supplementary figure paired with Figure 3 C/D: the number of genes in each
group's FOCUS 90% credible set.

    A  Genes in FOCUS 90% credible set per GTEx tissue   (39 tissues)
    B  Genes in FOCUS 90% credible set per OneK1K cell type (17 cell types)

Where C/D bin the credible-set genes by PIP for six representative groups, this
supplement gives the full count breakdown across every tissue / cell type:
one horizontal bar per group, sorted by descending credible-set size, with the
count printed at the bar end. Bars are coloured to match the A/B Miami panels
(royalblue tissues, cadetblue cell types).

Source data (each row is one credible-set gene-model; in_cred_set_pop1 == 1 for
all rows, so the per-group count = genes in that group's 90% credible set):
    ../data/focus_credset_gene_tissue_39.tsv   (tissues)
    ../data/focus_credset_gene_tissue_17CT.tsv  (cell types)

Run:  python figure_scripts/build_credset_count_panels.py
Out:  figures/panel_S_credset_tissue.png
      figures/panel_S_credset_celltype.png
      figures/figureS_credset_combined.png / .pdf
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

from build_pip_panels import pretty, pretty_ct  # shared label conventions

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("DATA_DIR", ROOT / "data"))
FIGS = ROOT / "figures"
plt.rcParams["font.family"] = "DejaVu Sans"

# match the A/B Miami significance colours
TISSUE_COLOR = "royalblue"
CT_COLOR     = "cadetblue"
FIG_H = 9.0          # shared height so the two panels align when placed side by side


def credset_counts(fname):
    """Genes in the FOCUS 90% credible set per group (rows == unique genes)."""
    df = pd.read_csv(DATA / fname, sep="\t")
    return df.groupby("tissues")["ens_gene_id"].nunique().sort_values(ascending=False)


def bar_panel(counts, letter, title, color, label_fn, save_path):
    fig, ax = plt.subplots(figsize=(7, FIG_H))
    fig.subplots_adjust(left=0.34, right=0.94, top=0.93, bottom=0.085)
    y = np.arange(len(counts))
    ax.barh(y, counts.values, color=color, edgecolor="white", linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels([label_fn(g) for g in counts.index], fontsize=10)
    ax.invert_yaxis()                       # largest credible set at the top
    ax.set_xlim(0, counts.max() * 1.12)
    ax.set_xlabel("Genes in FOCUS 90% credible set", fontsize=13)
    ax.tick_params(axis="x", labelsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    pad = counts.max() * 0.012
    for yi, v in zip(y, counts.values):
        ax.text(v + pad, yi, str(int(v)), ha="left", va="center", fontsize=9)

    ax.set_title(title, fontsize=14, fontweight="bold", loc="left", pad=12)
    fig.text(0.012, 0.985, letter, fontsize=30, fontweight="bold",
             ha="left", va="top")
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {save_path}  ({len(counts)} groups, "
          f"range {int(counts.min())}-{int(counts.max())})")


def _content_crop(img):
    a = np.asarray(img.convert("L"))
    ys, xs = np.where(a < 248)
    if len(xs) == 0:
        return img
    return img.crop((xs.min(), ys.min(), xs.max() + 1, ys.max() + 1))


def combine(left_png, right_png, out_png, out_pdf, gap=60, pad=24):
    """Place the two panels side by side (top-aligned) on a white canvas."""
    li = _content_crop(Image.open(left_png).convert("RGB"))
    ri = _content_crop(Image.open(right_png).convert("RGB"))
    h = max(li.height, ri.height)
    w = li.width + gap + ri.width + 2 * pad
    canvas = Image.new("RGB", (w, h + 2 * pad), (255, 255, 255))
    canvas.paste(li, (pad, pad))
    canvas.paste(ri, (pad + li.width + gap, pad))
    canvas.save(out_png)
    canvas.save(out_pdf, "PDF", resolution=300.0)
    print(f"saved {out_png}\n      {out_pdf}  ({canvas.width}x{canvas.height})")


def main():
    t = credset_counts("focus_credset_gene_tissue_39.tsv")
    c = credset_counts("focus_credset_gene_tissue_17CT.tsv")
    bar_panel(t, "A", "Credible-set genes per GTEx tissue",
              TISSUE_COLOR, pretty, FIGS / "panel_S_credset_tissue.png")
    bar_panel(c, "B", "Credible-set genes per OneK1K cell type",
              CT_COLOR, pretty_ct, FIGS / "panel_S_credset_celltype.png")
    combine(FIGS / "panel_S_credset_tissue.png",
            FIGS / "panel_S_credset_celltype.png",
            FIGS / "figureS_credset_combined.png",
            FIGS / "figureS_credset_combined.pdf")


if __name__ == "__main__":
    main()
