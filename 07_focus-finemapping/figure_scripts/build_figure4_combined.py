"""
Re-assemble combined Figure 4 from its four panel images.

Layout (panel letters read A->D top-to-bottom):
    A  Model construction and evaluation overview        (full width)
    B  Predicted-probability distributions on CAMP-Balanced   (full width)
    C  Pairwise AUC comparison        D  Top-quartile case enrichment   (split)

Panel sources (each is title-less; the bold letter + title are drawn here):
    A  panelA.png  -- the workflow schematic (six stages: variant-level ->
       single-feature gene-level -> cross-cohort features -> multi-tissue/
       cell-type unified PTRS -> multimodal -> CAMP evaluation), supplied
       directly as a tight-cropped PNG.
    B  figure4_panelB_violins.png
    C  figure4_panelC_pairwise.png
    D  figure4_panelD_or.png

Run:  python build_figure4_combined.py
Out:  figure4_combined.png / figure4_combined.pdf
"""
import os
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

FIGDIR = Path(os.environ.get("FIG_DIR", Path(__file__).resolve().parent))

PANELS = {
    "A": ("panelA.png",                  "Model construction and evaluation overview"),
    "B": ("figure4_panelB_violins.png",  "Predicted-probability distributions on CAMP-Balanced"),
    "C": ("figure4_panelC_pairwise.png", "Pairwise AUC comparison"),
    "D": ("figure4_panelD_or.png",       "Top-quartile case enrichment"),
}

# ---- layout in inches -------------------------------------------------------
M       = 0.20      # side margin
TOP     = 0.12      # top margin
BOT     = 0.12      # bottom margin
TB      = 0.50      # title band height (bold letter + title)
GAP     = 0.28      # vertical gap between panel blocks
CGAP    = 0.33      # horizontal gap between C and D
FIG_W   = 13.33
CONTENT_W = FIG_W - 2 * M

LETTER_FS = 19      # bold panel letter
TITLE_FS  = 14      # panel title
LETTER_DX = 0.30    # inches: title text x-offset past the letter


def aspect(path):
    im = mpimg.imread(path)
    h, w = im.shape[0], im.shape[1]
    return w / h, im


def main():
    imgs = {k: aspect(FIGDIR / fname) for k, (fname, _t) in PANELS.items()}

    aw, _ = imgs["A"]; bw, _ = imgs["B"]; cw, _ = imgs["C"]; dw, _ = imgs["D"]
    hA = CONTENT_W / aw
    hB = CONTENT_W / bw
    half = (CONTENT_W - CGAP) / 2.0
    hC = half / cw
    hD = half / dw
    hCD = max(hC, hD)

    FIG_H = TOP + TB + hA + GAP + TB + hB + GAP + TB + hCD + BOT
    fig = plt.figure(figsize=(FIG_W, FIG_H))

    def fx(v): return v / FIG_W
    def fy(v): return v / FIG_H

    def place(letter, left_in, bottom_in, w_in, h_in, title):
        ax = fig.add_axes([fx(left_in), fy(bottom_in), fx(w_in), fy(h_in)])
        ax.imshow(imgs[letter][1], aspect="auto")
        ax.axis("off")
        ty = bottom_in + h_in + 0.10          # baseline a touch above the panel
        fig.text(fx(left_in), fy(ty), letter, fontsize=LETTER_FS,
                 fontweight="bold", ha="left", va="bottom")
        fig.text(fx(left_in + LETTER_DX), fy(ty), title, fontsize=TITLE_FS,
                 ha="left", va="bottom")

    # rows from top
    yA = FIG_H - TOP - TB - hA
    place("A", M, yA, CONTENT_W, hA, PANELS["A"][1])
    yB = yA - GAP - TB - hB
    place("B", M, yB, CONTENT_W, hB, PANELS["B"][1])
    yCD = yB - GAP - TB - hCD
    place("C", M, yCD + (hCD - hC), half, hC, PANELS["C"][1])
    place("D", M + half + CGAP, yCD + (hCD - hD), half, hD, PANELS["D"][1])

    out_png = FIGDIR / "figure4_combined.png"
    out_pdf = FIGDIR / "figure4_combined.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)
    print(f"saved {out_png}\n      {out_pdf}  ({FIG_W:.2f}x{FIG_H:.2f} in)")


if __name__ == "__main__":
    main()
