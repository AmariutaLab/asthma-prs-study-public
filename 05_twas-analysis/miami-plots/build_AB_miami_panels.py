"""
Build ONLY the two Figure-3 Miami panels (A and B), not all 39 tissues / 17
cell types.

    Panel A : GTEx bulk-tissue TWAS Miami plot  -> tissue   `esophagus_mucosa`
    Panel B : OneK1K cell-type TWAS Miami plot  -> celltype `cd4_naive`

This is a trimmed copy of `miami_plot.ipynb`: same inputs, same Bonferroni
labelling threshold, same styling (royalblue tissue / cadetblue cell-type
highlights, Z=+/-1.96 dashed lines, adjustText repulsion). It just restricts
the per-context loop to the two contexts that become Fig-3 panels A and B, so
you can regenerate those without rendering the full set.

Outputs (SVG, dpi 1200; same filenames the notebook would produce so they stay
drop-in for the Affinity letter/title step):
    svg/miami_TISSUE_esophagus_mucosa.svg   (Panel A source)
    svg/miami_CELLTYPE_cd4_naive.svg        (Panel B source)

Run from this directory:
    python build_AB_miami_panels.py
"""
import os
from pathlib import Path
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from adjustText import adjust_text

# Editable text (not outlines) in SVG/PDF exports -- matches the notebook
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get('MIAMI_DATA_DIR', HERE / 'data'))
OUTDIR = Path(os.environ.get('MIAMI_SVG_DIR', HERE / 'svg'))
OUTDIR.mkdir(parents=True, exist_ok=True)

# The two contexts that become Fig-3 panels A / B
PANEL_A_TISSUE   = 'esophagus_mucosa'
PANEL_B_CELLTYPE = 'cd4_naive'


def generate_twas_miami_plot(
    df, name, significance_p, chrom_midpoints,
    output_fn="", significance_z=1.96,
    significance_color="royalblue", s=5, xaxis_padding=1.6e8,
    figsize=(15, 10), label_fontsize=14, axis_fontsize=18, tick_fontsize=14,
    context_fontsize=17,
):
    """Single-context Miami plot. Same data / colours / thresholds as
    miami_plot.ipynb, but with a larger canvas, enlarged fonts (label 14 / axis
    18 / tick 14 / context 17, up from 12.5/16/12/15 so the labels stay legible
    once the wide panel is shrunk into the Figure-3 grid, but small enough to
    de-overlap on the current canvas -- 16.5 was too big and left the chr17-22
    labels crowded), and a stronger, longer-running adjustText pass so the dense
    gene-label clusters (HLA on chr6, 17q21 on chr17) stop overlapping. adjustText
    re-fits the label boxes and positions to the bigger text; only the label
    LAYOUT differs from the notebook -- every plotted point, highlight and
    threshold line is identical."""
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    subset = df.loc[df["PANEL"] == name]
    name_formatted = subset["PANEL.CLEANED"].unique()[0]
    # Capitalize "naive"/"naïve" in the context label (e.g. CD4+ naïve T cells)
    name_formatted = name_formatted.replace("naïve", "Naïve").replace("naive", "Naive")
    highlighted = subset.loc[subset["TWAS.P"] < significance_p]

    ax.scatter(subset["Position"], subset["TWAS.Z"], c=subset["CHROM_COLOR"],
               s=s, alpha=.6, rasterized=True)
    ax.scatter(highlighted["Position"], highlighted["TWAS.Z"],
               c=significance_color, s=s + 2, alpha=.7, rasterized=True)
    ax.axhline(-significance_z, color="midnightblue", linestyle="--",
               label=f"Z = {significance_z}")
    ax.axhline(significance_z, color="midnightblue", linestyle="--",
               label=f"Z = {significance_z}")

    texts = []
    for _, row in highlighted.iterrows():
        texts.append(ax.text(
            row["Position"], row["TWAS.Z"], row["ID"],
            ha='center', va='center', fontsize=label_fontsize, fontstyle='italic',
            bbox=dict(facecolor='white', edgecolor=significance_color,
                      boxstyle='round,pad=.2', linewidth=.3)))

    # IMPORTANT: fix the axis limits BEFORE adjustText, with a horizontal gutter
    # (xaxis_padding) and vertical headroom (zpad). ensure_inside_axes then keeps
    # every label inside this padded box, so labels settle in the gutter instead
    # of flush against the y-axis / top / bottom spines.
    zmin, zmax = subset["TWAS.Z"].min(), subset["TWAS.Z"].max()
    zspan = zmax - zmin
    ax.set_xlim(df["Position"].min() - xaxis_padding,
                df["Position"].max() + xaxis_padding)
    # Extra top headroom so (a) the context label sits flush under the top box
    # line with a clear strip and (b) the top gene labels are not crammed.
    ax.set_ylim(zmin - 0.12 * zspan, zmax + 0.18 * zspan)

    # Context label (tissue / cell-type name) flush against the upper box line.
    ctx_text = ax.text(0.012, 0.985, name_formatted, transform=ax.transAxes,
                       ha='left', va='top', fontsize=context_fontsize)

    # Stronger, longer separation than the notebook: force_explode first blows
    # apart any initially-overlapping boxes, larger expand/force_text keep them
    # apart, and the raised iter/time limits let the solver fully converge on
    # the dense HLA / 17q21 clusters. ensure_inside_axes keeps labels in-box.
    adjust_text(
        texts,
        objects=[ctx_text],          # keep gene labels off the context label
        expand=(2.3, 2.7),
        force_text=(0.95, 1.3),
        force_explode=(0.78, 1.0),
        force_static=(0.12, 0.2),
        max_move=550,
        ensure_inside_axes=True,
        iter_lim=8000,
        time_lim=220,
        arrowprops=dict(arrowstyle='-', color='#A3C5D9', lw=0.6),
        min_arrow_len=2,
        ax=ax,
    )

    ax.set_xticks(chrom_midpoints)
    ax.set_xticklabels(chrom_midpoints.index.astype(int), fontsize=tick_fontsize)
    ax.tick_params(axis='y', labelsize=tick_fontsize)
    ax.set_xlabel("Chromosome", fontsize=axis_fontsize)
    ax.set_ylabel("TWAS Z-score", fontsize=axis_fontsize)

    plt.tight_layout()
    if output_fn:
        plt.savefig(output_fn, format="svg", dpi=1200)
        print(f"  wrote {output_fn}")
    plt.close(fig)


def prep(df):
    """Add Position / chromosome-colour columns exactly as the notebook does."""
    chrom_padding = 4e7
    chrom_lengths = df.sort_values(["CHR", "P1"]).groupby("CHR")["P1"].max() + chrom_padding
    chrom_lengths[22] = chrom_lengths[22] - chrom_padding  # last chr needs no pad
    chrom_offsets = chrom_lengths.cumsum().shift(1, fill_value=0)
    df["Position"] = df["P1"] + df["CHR"].map(chrom_offsets)
    midpoints = df.groupby("CHR")["Position"].mean()
    chrom_colors = ["silver", "dimgray"]
    cmap = {c: chrom_colors[i % 2] for i, c in enumerate(df["CHR"].unique())}
    df["CHROM_COLOR"] = df["CHR"].map(cmap)
    return df, midpoints


def main():
    df_t = pd.read_csv(DATA / "twas_all-tissue-all-chr.tsv", sep="\t",
                       index_col=0).dropna(subset=["TWAS.Z", "TWAS.P"])
    df_c = pd.read_csv(DATA / "twas_all-celltype-all-chr.tsv", sep="\t",
                       index_col=0).dropna(subset=["TWAS.Z", "TWAS.P"])

    # Bonferroni thresholds over the full (pre-QC) test universe -- unchanged
    bonf_t = 0.05 / df_t.shape[0]
    bonf_c = 0.05 / df_c.shape[0]
    print(f"Tissue   pairs={df_t.shape[0]:,}  Bonferroni p<{bonf_t:.3e}")
    print(f"Celltype pairs={df_c.shape[0]:,}  Bonferroni p<{bonf_c:.3e}")

    df_t, mid_t = prep(df_t)
    df_c, mid_c = prep(df_c)

    # Panel A -- tissue (royalblue highlights). Esophagus mucosa has ~90
    # significant genes, so it gets a larger canvas to fit the bigger font
    # without overlaps.
    print(f"Panel A: tissue {PANEL_A_TISSUE}")
    generate_twas_miami_plot(
        df_t, PANEL_A_TISSUE, bonf_t, mid_t, figsize=(21, 13),
        output_fn=str(OUTDIR / f"miami_TISSUE_{PANEL_A_TISSUE}.svg"))

    # Panel B -- cell type (cadetblue highlights, matching the notebook)
    print(f"Panel B: celltype {PANEL_B_CELLTYPE}")
    generate_twas_miami_plot(
        df_c, PANEL_B_CELLTYPE, bonf_c, mid_c,
        significance_color="cadetblue",
        output_fn=str(OUTDIR / f"miami_CELLTYPE_{PANEL_B_CELLTYPE}.svg"))


if __name__ == "__main__":
    main()
