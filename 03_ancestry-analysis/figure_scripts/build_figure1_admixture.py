"""
Figure 1, Panel B -- ADMIXTURE ancestry-proportion plot.

Seven side-by-side stacked-bar panels (one bar per individual) of the K=5
ADMIXTURE proportions (AFR, AMR, EAS, EUR, SAS): the study cohorts GACRS and
CAMP, followed by the five 1000 Genomes super-populations as references. Within
each panel individuals are sorted by descending proportion, and panel widths are
scaled to each cohort's sample count.

Input (see CONFIG):
    <ADMIXTURE_CSV> -- long/wide per-individual proportions with columns:
        IID, Population, "1000G AFR", "1000G AMR", "1000G EAS",
        "1000G EUR", "1000G SAS"
      Population values: GACRS, CAMP, 1000G AFR, 1000G AMR, 1000G EAS,
        1000G EUR, 1000G SAS.

    !! INDIVIDUAL-LEVEL DATA !!  The CSV is keyed by individual sample ID (IID)
    for the GACRS and CAMP cohorts, so it is NOT committed to the public
    repository -- it lives only in the private repo. This script (which merely
    reads it) is safe to release; the input is not.

Output (see CONFIG -> FIG_OUTDIR, default = ../figures):
    figure1_admixture.svg   (editable text, dpi 1200)
    figure1_admixture.png   (300 dpi preview)
"""
import os
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"


# =============================================================================
# CONFIG
# =============================================================================
DATA_ROOT = Path(os.environ.get(
    "FIG1_DATA_ROOT",
    "/Users/nancyh/Desktop/hartwell-asthma-main/local-computer-files/figure-1/data"))

ADMIXTURE_CSV = Path(os.environ.get(
    "ADMIXTURE_CSV", DATA_ROOT / "admixture-results" /
    "02_gacrs-camp-1kg-cohort-admixture-visualization_admixture-proportion-all-pop.csv"))

FIG_OUTDIR = Path(os.environ.get(
    "FIG_OUTDIR", Path(__file__).resolve().parents[1] / "figures"))
FIG_OUTDIR.mkdir(parents=True, exist_ok=True)

# AFR, AMR, EAS, EUR, SAS  (yellow, red, green, royal-blue, indigo)
PALETTE = ["#f9cc29", "#d24040", "#39853d", "#4169E1", "#241882"]
DATA_COLUMNS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
# Panels left-to-right: study cohorts first, then the 1000G reference pops
COHORTS = ["GACRS", "CAMP", "1000G AFR", "1000G AMR", "1000G EAS", "1000G EUR", "1000G SAS"]
XLABELS = ["GACRS", "CAMP", "AFR", "AMR", "EAS", "EUR", "SAS"]
BAR_WIDTH = 2


# =============================================================================
# Load + split per cohort
# =============================================================================
if not ADMIXTURE_CSV.exists():
    raise FileNotFoundError(
        f"ADMIXTURE proportions CSV not found: {ADMIXTURE_CSV}\n"
        f"This is individual-level data (private repo only); set ADMIXTURE_CSV "
        f"to point at your copy.")

df = pd.read_csv(ADMIXTURE_CSV)

dfs = []
for pop in COHORTS:
    sub = df.loc[df["Population"] == pop].rename(
        columns={c: c.replace("1000G ", "") for c in df.columns})
    sub = sub.set_index("IID")[DATA_COLUMNS].sort_values(DATA_COLUMNS, ascending=False)
    dfs.append(sub)
    print(f"  {pop:10s} : {sub.shape[0]:,} individuals")


# =============================================================================
# Render -- 7 panels, widths proportional to sample counts
# =============================================================================
widths = {"width_ratios": [sub.shape[0] for sub in dfs]}
fig, axes = plt.subplots(1, 7, figsize=(9.25, 3), sharey=True, gridspec_kw=widths)

for sub, ax, xlabel in zip(dfs, axes.flatten(), XLABELS):
    g = sub.plot.bar(stacked=True, width=BAR_WIDTH, edgecolor=None, linewidth=1,
                     color=PALETTE, ax=ax, rasterized=True)
    g.set_xlabel(xlabel)
    g.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    if xlabel == "GACRS":
        g.set_ylim(0, 1)
        g.set_ylabel("ADMIXTURE Proportion")
    else:
        g.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    if xlabel != "SAS":
        g.get_legend().set_visible(False)
    else:
        g.legend(ncols=1, loc="lower left", bbox_to_anchor=(1, .5),
                 labelspacing=0.5, frameon=False)

plt.subplots_adjust(wspace=.05)

out_svg = FIG_OUTDIR / "figure1_admixture.svg"
out_png = FIG_OUTDIR / "figure1_admixture.png"
plt.savefig(out_svg, format="svg", dpi=1200, bbox_inches="tight")
plt.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved: {out_svg}\n       {out_png}")
