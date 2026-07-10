"""
Render the committed Miami SVGs into the Figure-3 panel A / B PNGs.

The Miami panels are produced upstream as SVGs by
`05_twas-analysis/miami-plots/build_AB_miami_panels.py` (de-overlapped labels,
"TWAS Z-score" axis, royalblue/cadetblue highlights). This script turns those
SVGs into the letter+title-bearing raster panels that build_figure3_combined.py
places:

    svg/miami_TISSUE_esophagus_mucosa.svg -> figures/panel_A_miami_tissue.png   (A)
    svg/miami_CELLTYPE_cd4_naive.svg      -> figures/panel_B_miami_celltype.png  (B)

Pipeline per panel:
  1. Expand the SVG's viewBox/width/height by a large margin. adjustText flings
     some labels (e.g. the chr17 GSDMB/ORMDL3 cluster) well past the figure edge;
     matplotlib's tight_layout leaves no room, so a renderer clips them. The
     margin gives them room; content_crop trims the excess whitespace back.
  2. Render to PNG with headless Chrome (the only SVG renderer available here;
     cairosvg/rsvg/qlmanage were not usable). Scale factor 1 keeps the screenshot
     small enough that Chrome does not OOM-crash on the padded canvas.
  3. Content-crop to the real plot, then compose under a bold panel letter +
     title band (DejaVu Sans Bold, matching panels C-F).

Requires Google Chrome at the standard macOS path. Run:
    python figure_scripts/render_miami_AB_from_svg.py
"""
import os
import re
import shutil
import tempfile
import subprocess
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SVGDIR = Path(os.environ.get(
    "MIAMI_SVG_DIR", ROOT.parent / "05_twas-analysis" / "miami-plots" / "svg"))
FIGDIR = ROOT / "figures"
CHROME = os.environ.get(
    "CHROME_BIN", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
plt.rcParams["font.family"] = "DejaVu Sans"

# (svg name, (svg width pt, height pt), letter, title, output png)
PANELS = [
    ("miami_TISSUE_esophagus_mucosa.svg", "A",
     "GTEx bulk-tissue TWAS Miami plot", "panel_A_miami_tissue.png"),
    ("miami_CELLTYPE_cd4_naive.svg", "B",
     "OneK1K cell-type TWAS Miami plot", "panel_B_miami_celltype.png"),
]

MARGIN     = 1200   # SVG units added on every side so flung labels are not clipped
SCALE      = 1      # Chrome device-scale-factor (2/3 OOM-crash on the padded canvas)
FONT_SCALE = 1.0    # The committed SVGs are now rebuilt by build_AB_miami_panels.py
                    # with the larger fonts already baked in (label 16.5 / axis 21 /
                    # tick 16 / context 20) and adjustText-fitted boxes, so no in-place
                    # scaling is needed. (Was 1.32 as a stopgap for the old 12.5px SVGs.)
LETTER_FS  = 40
TITLE_FS   = 22
FIG_W      = 12.0   # inches
TITLE_BAND = 0.66   # inches reserved on top for the letter + title


def scale_fonts(txt, scale=FONT_SCALE):
    """Multiply every `font-size: Npx` in the SVG by `scale` so the whole panel's
    text is larger. Boxes/positions are left as adjustText placed them."""
    if scale == 1.0:
        return txt
    return re.sub(r"font-size:\s*([\d.]+)px",
                  lambda m: f"font-size: {float(m.group(1)) * scale:.3g}px", txt)


def expand_svg(svg, out_svg, margin=MARGIN):
    """Scale up fonts, then grow the <svg> box by `margin` on every side so flung
    labels are not clipped; return the new (w, h) pt."""
    txt = scale_fonts(svg.read_text())
    tag = re.search(r"<svg[^>]*>", txt).group(0)
    w = float(re.search(r'width="([\d.]+)pt"', tag).group(1))
    h = float(re.search(r'height="([\d.]+)pt"', tag).group(1))
    x0, y0, vw, vh = map(float, re.search(
        r'viewBox="([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+)"', tag).groups())
    new = re.sub(r'width="[\d.]+pt"',  f'width="{w + 2*margin:g}pt"', tag)
    new = re.sub(r'height="[\d.]+pt"', f'height="{h + 2*margin:g}pt"', new)
    new = re.sub(r'viewBox="[-\d. ]+"',
                 f'viewBox="{x0-margin:g} {y0-margin:g} {vw+2*margin:g} {vh+2*margin:g}"',
                 new)
    out_svg.write_text(txt.replace(tag, new, 1))
    return int(round(w + 2*margin)), int(round(h + 2*margin))


def render(svg, raw_png, tmp):
    padded = tmp / (svg.stem + ".padded.svg")
    win = expand_svg(svg, padded)
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-sandbox",
                    "--hide-scrollbars", f"--force-device-scale-factor={SCALE}",
                    f"--window-size={win[0]},{win[1]}",
                    "--default-background-color=FFFFFFFF",
                    f"--screenshot={raw_png}", f"file://{padded}"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def content_crop(path, thresh=250):
    im = Image.open(path).convert("RGB")
    a = np.asarray(im.convert("L"))
    ys, xs = np.where(a < thresh)
    return im.crop((xs.min(), ys.min(), xs.max() + 1, ys.max() + 1))


def compose(raw_png, letter, title, out_png):
    plot = content_crop(raw_png)
    pw, ph = plot.size
    plot_h_in = FIG_W * ph / pw
    fig_h = plot_h_in + TITLE_BAND
    fig = plt.figure(figsize=(FIG_W, fig_h))
    ax = fig.add_axes([0.0, 0.0, 1.0, plot_h_in / fig_h])
    ax.imshow(plot); ax.axis("off")
    ty = 1.0 - 0.12 * TITLE_BAND / fig_h
    fig.text(0.004, ty, letter, fontsize=LETTER_FS, fontweight="bold",
             ha="left", va="top")
    fig.text(0.052, ty - 0.004, title, fontsize=TITLE_FS, fontweight="bold",
             ha="left", va="top")
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"  {out_png.name}: content {pw}x{ph} -> fig {FIG_W}x{fig_h:.2f}in")


def main():
    if not Path(CHROME).exists():
        raise SystemExit(f"Chrome not found at {CHROME} (set CHROME_BIN)")
    tmp = Path(tempfile.mkdtemp(prefix="miami_AB_"))
    try:
        for svgname, letter, title, outname in PANELS:
            svg = SVGDIR / svgname
            if not svg.exists():
                raise SystemExit(f"missing SVG: {svg}")
            raw = tmp / f"{letter}_raw.png"
            render(svg, raw, tmp)
            compose(raw, letter, title, FIGDIR / outname)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("done")


if __name__ == "__main__":
    main()
