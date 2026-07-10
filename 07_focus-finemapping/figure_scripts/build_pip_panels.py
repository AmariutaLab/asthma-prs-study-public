"""
Build Figure 3 panels C and D: per-group stacked-proportion bars of FOCUS PIP.

    C  Gene-tissue pairs by FOCUS PIP (GTEx tissues)
    D  Gene-cell-type pairs by FOCUS PIP (OneK1K cell types)

Each panel: one horizontal stacked bar per group (tissue / cell type), rows
ordered by descending count of PIP > 0.1, with a "PIP > 0.1 / N=" annotation in
a fixed right-margin column and the PIP-bin legend above the bars. A big bold
panel letter sits at the top-left and the title is padded clear of the bars.

Binning (corrected): right-CLOSED bins so the boundary at 0.1 is consistent with
the "PIP > 0.1" count -- a PIP of exactly 0.1 falls in the 0-0.1 (low) bin, not
the 0.1-0.3 segment. (The old notebook used right=False, which put PIP==0.1 in
the colored 0.1-0.3 segment while excluding it from the N>0.1 count.)

Source data: ../data/focus_credset_gene_tissue_39.tsv (tissues),
             ../data/focus_credset_gene_tissue_17CT.tsv (cell types).
Run:  python figure_scripts/build_pip_panels.py
Out:  figures/panel_C_pip_tissue.png, figures/panel_D_pip_celltype.png
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("DATA_DIR", ROOT / "data"))
FIGS = ROOT / "figures"
plt.rcParams["font.family"] = "DejaVu Sans"

# right-CLOSED edges; with include_lowest the first bin is [0,0.1].
PIP_BINS   = [0.0, 0.1, 0.3, 0.5, 0.8, 1.0]
PIP_LABELS = ["0-0.1", "0.1-0.3", "0.3-0.5", "0.5-0.8", "0.8-1.0"]
PIP_COLORS = ["#F0E442", "#56B4E9", "#009E73", "#E69F00", "#D55E00"]
PIP_HIGH   = 0.1

SELECTED_TISSUES = ["thyroid", "esophagus_mucosa", "heart_atrial_appendage",
                    "adipose_subcutaneous", "colon_transverse", "artery_coronary"]
SELECTED_CTS     = ["nk", "cd4_tcm", "cd4_naive", "nk_cd56bright", "cd4_ctl", "treg"]

# Cell-type display names follow the standard immunology / OneK1K convention
# (NK, CD4, CD8, MAIT as fixed acronyms; memory/effector subsets as TCM/TEM/TCM;
# subset modifiers like "Naive"/"CD56bright"/"gdT" in their conventional case)
# rather than naive Title-casing ("Nk", "Cd4 Tcm", "Gdt").
CT_LABELS = {
    "nk":             "NK",
    "nk_cd56bright":  "NK CD56bright",
    "cd4_naive":      "CD4 Naive",
    "cd4_tcm":        "CD4 TCM",
    "cd4_tem":        "CD4 TEM",
    "cd4_ctl":        "CD4 CTL",
    "cd8_naive":      "CD8 Naive",
    "cd8_tcm":        "CD8 TCM",
    "cd8_tem":        "CD8 TEM",
    "treg":           "Treg",
    "b_naive":        "B Naive",
    "b_memory":       "B Memory",
    "b_intermediate": "B Intermediate",
    "cd14_mono":      "CD14 Mono",
    "cd16_mono":      "CD16 Mono",
    "mait":           "MAIT",
    "gdt":            "gdT",
}


def pretty(name):
    return name.replace("_", " ").title()


def pretty_ct(name):
    """Cell-type labels in standard immunology convention (see CT_LABELS)."""
    return CT_LABELS.get(name, pretty(name))


def load(fname):
    df = pd.read_csv(DATA / fname, sep="\t")
    df["pips_pop1"] = pd.to_numeric(df["pips_pop1"], errors="coerce")
    df.dropna(subset=["pips_pop1", "tissues", "symbol"], inplace=True)
    return df


def pip_panel(df, selected, letter, title, save_path, label_fn=pretty):
    sub = df[df["tissues"].isin(selected)].copy()

    # order rows + annotate by count of PIP > 0.1 (strict)
    counts = (sub[sub["pips_pop1"] > PIP_HIGH]["tissues"].value_counts()
              .reindex(selected, fill_value=0).sort_values(ascending=False))
    order = counts.index.tolist()

    # right-closed bins -> PIP==0.1 lands in 0-0.1, consistent with the N>0.1 count
    sub["pip_bin"] = pd.cut(sub["pips_pop1"], bins=PIP_BINS, labels=PIP_LABELS,
                            include_lowest=True, right=True)
    bin_counts = (sub.groupby(["tissues", "pip_bin"], observed=False).size()
                  .unstack(fill_value=0)
                  .reindex(index=order, columns=PIP_LABELS, fill_value=0))
    props = bin_counts.div(bin_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    # Proportions matched to the original styled panel (title spans ~75% width,
    # big bold letter at the corner). A reserved top band holds the letter +
    # title + legend, clear of the bars; tight bbox keeps long labels from
    # clipping.
    fig, ax = plt.subplots(figsize=(9, 4.8))
    fig.subplots_adjust(left=0.205, right=0.865, top=0.72, bottom=0.135)
    y, left = np.arange(len(order)), np.zeros(len(order))
    for label, color in zip(PIP_LABELS, PIP_COLORS):
        ax.barh(y, props[label].values, left=left, color=color,
                label=f"PIP {label}", edgecolor="white", linewidth=0.7)
        left += props[label].values

    for yi, cnt in zip(y, counts.values):
        ax.text(1.025, yi, f"PIP > {PIP_HIGH}\nN={int(cnt)}", ha="left", va="center",
                fontsize=10, transform=ax.get_yaxis_transform())

    ax.set_yticks(y)
    ax.set_yticklabels([label_fn(t) for t in order], fontsize=13)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xticks(np.arange(0, 1.01, 0.25))
    ax.set_xlabel("Proportion of gene pairs", fontsize=14)
    ax.tick_params(axis="x", labelsize=11)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # PIP-bin legend: horizontal strip just above the bars, below the title.
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.04), ncol=5, fontsize=10,
              frameon=False, handlelength=1.1, columnspacing=1.1,
              handletextpad=0.4, borderaxespad=0.0)

    # Big bold panel letter at the top-left corner; title on the same line to its
    # right (left-aligned), both clear of the legend/bars below.
    fig.text(-0.043, 0.965, letter, fontsize=33, fontweight="bold",
             ha="left", va="top")
    fig.text(0.005, 0.935, title, fontsize=18, fontweight="bold",
             ha="left", va="top")

    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {save_path}  (N>0.1 per row: {dict(zip(order, counts.values))})")
    return counts


def main():
    tissue_df = load("focus_credset_gene_tissue_39.tsv")
    ct_df = load("focus_credset_gene_tissue_17CT.tsv")
    pip_panel(tissue_df, SELECTED_TISSUES, "C",
              "Credible-set gene-tissue pairs by FOCUS PIP (GTEx tissues)",
              FIGS / "panel_C_pip_tissue.png")
    pip_panel(ct_df, SELECTED_CTS, "D",
              "Credible-set gene-cell-type pairs by FOCUS PIP (OneK1K cell types)",
              FIGS / "panel_D_pip_celltype.png", label_fn=pretty_ct)


if __name__ == "__main__":
    main()
