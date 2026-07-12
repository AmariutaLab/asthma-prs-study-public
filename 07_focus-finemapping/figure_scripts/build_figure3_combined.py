"""
Assemble combined Figure 3 (panels A-F) into a 3x2 grid.

    A  GTEx bulk-tissue TWAS Miami plot        B  OneK1K cell-type TWAS Miami plot
    C  Gene-tissue pairs by FOCUS PIP          D  Gene-cell-type pairs by FOCUS PIP
    E  Regional FOCUS: 17q21 / GSDMB (GTEx)    F  Regional FOCUS: HINT1 (OneK1K)

Panel content is laid out so the letters read A->F top-to-bottom by row: the
PIP-bar panels are C/D (middle row) and the regional FOCUS panels are E/F
(bottom row). Each panel PNG carries its own baked letter + title, so this
script only places them.

Panel sources / how to refresh a letter-bearing input:
  A/B  Miami plots -- styled, letter-bearing fixed inputs from the original
       Affinity assembly: `panel_A_miami_tissue.png`, `panel_B_miami_celltype.png`.
  C/D  PIP-bar panels -- regenerated from data by build_pip_panels.py
       (`panel_C_pip_tissue.png`, `panel_D_pip_celltype.png`); right-closed PIP
       bins, big corner letter + title.
  E/F  Regional FOCUS panels -- regenerate from regional_focus_plot.py with
       `--panel-letter E`/`F`: `panel_E_17q21_tissues.png`, `panel_F_celltypes.png`.
       (The repo/README also keep the C/D-lettered standalone versions for
       single-panel use; only the combined figure uses the E/F-lettered ones.)

Layout boxes were measured from the published figure3_combined_ABCDEF.png by
whitespace-gap detection, so the output matches that geometry.
"""
import os
from pathlib import Path
from PIL import Image

FIGDIR = Path(os.environ.get(
    "FIG_DIR", Path(__file__).resolve().parents[1] / "figures"))

CANVAS = (2471, 2418)            # (width, height) of the combined figure
# (left, top, right, bottom) cell box for each panel, in canvas pixels.
# Row heights match each row's panel aspect so there is no empty band: the A/B
# Miami panels are wider-than-tall (the de-overlapped 21x13 / 15x10 renders),
# so their row is shorter than the older near-square version. Inter-row gaps are
# kept tight (~36 px) so the wide A/B panels get as much height as possible.
PANELS = [
    ("panel_A_miami_tissue.png",      (24, 18, 1218, 880)),
    ("panel_B_miami_celltype.png",    (1270, 18, 2452, 880)),
    # C/D (PIP-bar panels) occupy the middle row. The panels carry a white
    # spacer below the title (added in the relabel step) so the title sits clear
    # of the legend/bars; the row height matches their padded aspect.
    ("panel_C_pip_tissue.png",        (24, 916, 1218, 1516)),
    ("panel_D_pip_celltype.png",      (1270, 916, 2452, 1516)),
    # ... and E/F (regional FOCUS panels) the bottom row, sized to fit snugly
    # so they sit close under C/D.
    ("panel_E_17q21_tissues.png",     (24, 1552, 1218, 2402)),
    ("panel_F_celltypes.png",         (1270, 1552, 2452, 2402)),
]
OUT_PNG = FIGDIR / "figure3_combined_ABCDEF.png"
OUT_PDF = FIGDIR / "figure3_combined_ABCDEF.pdf"
PDF_DPI = 300


def _content_crop(img):
    """Crop away the panel's surrounding whitespace so its real content (letter,
    title, plot) starts at the image's top-left -- this normalizes the different
    internal margins of the Affinity (A/B) vs matplotlib (C-F) panels."""
    import numpy as np
    a = np.asarray(img.convert("L"))
    ys, xs = np.where(a < 248)
    if len(xs) == 0:
        return img
    return img.crop((xs.min(), ys.min(), xs.max() + 1, ys.max() + 1))


def _fit_height(img, box, pad=6):
    """Height the (already content-cropped) image would take if fit into its box
    preserving aspect (min of width- and height-limited scaling)."""
    l, t, r, b = box
    cw, ch = (r - l) - 2 * pad, (b - t) - 2 * pad
    return img.height * min(cw / img.width, ch / img.height)


def place(canvas, img, box, target_h, pad=6):
    """Scale the image to `target_h` (preserving aspect) and paste at the box's
    TOP-LEFT (+pad). Panels are processed in left/right row pairs and given the
    SAME target height, so the two panels in a row share top and bottom edges
    (their letters stay left-aligned per column)."""
    l, t, r, b = box
    s = target_h / img.height
    nw, nh = max(1, int(round(img.width * s))), max(1, int(round(img.height * s)))
    img2 = img.resize((nw, nh), Image.LANCZOS)
    canvas.paste((255, 255, 255), (l, t, r, b))
    canvas.paste(img2, (l + pad, t + pad))
    return nw, nh


def main():
    canvas = Image.new("RGB", CANVAS, (255, 255, 255))
    imgs = {}
    for name, _box in PANELS:
        path = FIGDIR / name
        if not path.exists():
            raise SystemExit(f"missing panel image: {path}")
        imgs[name] = _content_crop(Image.open(path).convert("RGB"))
    # Process in left/right row pairs; align each pair to the SMALLER natural fit
    # height so the panel already filling its column width (e.g. A, D) is unchanged
    # and its taller-aspect neighbour is shrunk to match -> bottoms line up.
    for i in range(0, len(PANELS), 2):
        (nl, bl), (nr, br) = PANELS[i], PANELS[i + 1]
        target_h = min(_fit_height(imgs[nl], bl), _fit_height(imgs[nr], br))
        for name, box in (PANELS[i], PANELS[i + 1]):
            nw, nh = place(canvas, imgs[name], box, target_h)
            print(f"  {name:34s} -> box {box}  pasted {nw}x{nh}")
    canvas.save(OUT_PNG)
    canvas.save(OUT_PDF, "PDF", resolution=float(PDF_DPI))
    print(f"saved {OUT_PNG}\n      {OUT_PDF}  ({CANVAS[0]}x{CANVAS[1]})")


if __name__ == "__main__":
    main()
