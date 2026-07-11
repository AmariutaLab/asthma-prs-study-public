"""
Build a single-page combined Figure4.pdf:

  Row 1  Panel A  (from the user's slide 1 capture — preserves BioRender icons)
  Row 2  Panel B  (matplotlib violins)
  Row 3  Panel C + Panel D  (side-by-side)

Layout rules:
  * Tight inter-row padding (no big white gaps)
  * Panel labels "A   Title", "B   …", "C   …", "D   …" at top-left of each
    panel, A/B/C aligned at the SAME left x, D at the start of its right column
  * Helvetica bold for labels, 14 pt at the final 16-inch wide canvas (scales to
    ~6 pt at Cell 6.85-inch print width — author can scale up later if needed)
  * Output: Figure4.pdf, 300 dpi
"""
import os
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont, ImageOps

FIGDIR = Path(os.environ.get('FIG4_ROOT', '/Users/nancyh/Desktop/hartwell/gene_model/score/combine')) / 'figures'  # override with FIG4_ROOT env var
OUT_PDF = FIGDIR / 'Figure4.pdf'

# --- Sources ----------------------------------------------------------------
# Prefer a user-provided high-quality Panel A export from PowerPoint if it exists,
# otherwise fall back to the macOS Quick Look thumbnail (which may show z-order
# / text-overlap artifacts that aren't actually in your slide).
PANEL_A_CANDIDATES = [
    FIGDIR / 'panelA.png',                         # curated Panel A — the only source we use
]
PANEL_A_PNG = None
for cand in PANEL_A_CANDIDATES:
    if cand.exists():
        PANEL_A_PNG = cand
        print(f"Using Panel A from: {PANEL_A_PNG}")
        from PIL import Image as _IM
        _w, _h = _IM.open(PANEL_A_PNG).size
        print(f"  Source dimensions: {_w} × {_h}  (target canvas ~ 4000 px wide)")
        if _w < 3000:
            print(f"  NOTE: source is lower-resolution than ideal for print. For best Cell")
            print(f"        print quality (300 dpi @ 6.85 in = 2055 px wide), re-export from")
            print(f"        PowerPoint using Preferences → Save → Default resolution = 300 dpi.")
        break
if PANEL_A_PNG is None:
    raise FileNotFoundError("No Panel A image found. Export Slide 1 from PowerPoint.")
PANEL_B_PNG = FIGDIR / 'figure4_panelB_violins.png'
PANEL_C_PNG = FIGDIR / 'figure4_panelC_pairwise.png'
PANEL_D_PNG = FIGDIR / 'figure4_panelD_or.png'

# --- Layout knobs (px, at 256 dpi target = 4096 px wide canvas ~ 16 in) -----
W = 4096                  # canvas width
DPI = 256                 # 4096 / 16
PAD_HORIZ = 30            # horizontal pad between panels
PAD_VERT  = 12            # vertical pad between rows (tightened)
LABEL_HEIGHT = 110        # space reserved for each panel label
LABEL_LEFT_OFFSET = 24    # x-offset of label from panel-area left edge
LABEL_FONT_SZ = 80        # ~28 pt at 256 dpi — Cell allows up to ~12 pt at
                          # 6.85-in print → 80 px scales down to ~13 pt
LABEL_FONT_PATH_CANDIDATES = [
    '/System/Library/Fonts/Helvetica.ttc',
    '/Library/Fonts/Arial.ttf',
    '/System/Library/Fonts/Supplemental/Arial.ttf',
]

label_font = None
for p in LABEL_FONT_PATH_CANDIDATES:
    if Path(p).exists():
        try:
            label_font = ImageFont.truetype(p, LABEL_FONT_SZ)
            print(f"Using font: {p}")
            break
        except Exception:
            continue
if label_font is None:
    label_font = ImageFont.load_default()
    print("WARNING: falling back to PIL default font")


# --- Load and prep Panel A ---------------------------------------------------
panel_a = Image.open(PANEL_A_PNG).convert('RGB')
print(f"Panel A original: {panel_a.size}")

# Auto-trim white border so Panel A starts at content (max ~5% off any side)
def trim_white(img, threshold=250, max_pct=0.10):
    """Crop surrounding white margins (within max_pct of dimensions)."""
    arr = img.load()
    w, h = img.size
    mx, my = int(w * max_pct), int(h * max_pct)
    # find content bbox by scanning rows/cols
    def is_white_row(y):
        return all(sum(arr[x, y]) > threshold * 3 - 10 for x in range(0, w, 50))
    def is_white_col(x):
        return all(sum(arr[x, y]) > threshold * 3 - 10 for y in range(0, h, 50))
    top = 0
    while top < my and is_white_row(top): top += 1
    bottom = h - 1
    while bottom > h - my and is_white_row(bottom): bottom -= 1
    left = 0
    while left < mx and is_white_col(left): left += 1
    right = w - 1
    while right > w - mx and is_white_col(right): right -= 1
    return img.crop((left, max(0, top - 10), right + 1, min(h, bottom + 10)))

panel_a = trim_white(panel_a, max_pct=0.75)
print(f"Panel A after trim:  {panel_a.size}")

# Resize Panel A to canvas width minus side padding
new_w_a = W - 2 * PAD_HORIZ
ratio_a = new_w_a / panel_a.width
panel_a = panel_a.resize((new_w_a, int(panel_a.height * ratio_a)),
                         Image.LANCZOS)
print(f"Panel A resized:     {panel_a.size}")

# --- Load Panels B/C/D and resize -------------------------------------------
panel_b = Image.open(PANEL_B_PNG).convert('RGB')
new_w_b = W - 2 * PAD_HORIZ
panel_b = panel_b.resize((new_w_b, int(panel_b.height * new_w_b / panel_b.width)),
                         Image.LANCZOS)
print(f"Panel B resized:     {panel_b.size}")

half_w = (W - 2 * PAD_HORIZ - PAD_HORIZ) // 2  # subtract gap between C and D
panel_c = Image.open(PANEL_C_PNG).convert('RGB')
panel_c = panel_c.resize((half_w, int(panel_c.height * half_w / panel_c.width)),
                         Image.LANCZOS)
panel_d = Image.open(PANEL_D_PNG).convert('RGB')
panel_d = panel_d.resize((half_w, int(panel_d.height * half_w / panel_d.width)),
                         Image.LANCZOS)
# normalize C and D to the same height
cd_h = max(panel_c.height, panel_d.height)
def pad_to_height(img, h):
    if img.height == h: return img
    new = Image.new('RGB', (img.width, h), 'white')
    new.paste(img, (0, (h - img.height) // 2))
    return new
panel_c = pad_to_height(panel_c, cd_h)
panel_d = pad_to_height(panel_d, cd_h)
print(f"Panel C/D resized:   {panel_c.size} / {panel_d.size}")


# --- Compose canvas ----------------------------------------------------------
total_h = (
    PAD_VERT                                      # top
    + LABEL_HEIGHT + panel_a.height               # A row
    + PAD_VERT
    + LABEL_HEIGHT + panel_b.height               # B row
    + PAD_VERT
    + LABEL_HEIGHT + cd_h                         # C+D row
    + PAD_VERT                                    # bottom
)
canvas = Image.new('RGB', (W, total_h), 'white')
draw = ImageDraw.Draw(canvas)

LABEL_X_LEFT = PAD_HORIZ + LABEL_LEFT_OFFSET  # consistent x for A, B, C
LABEL_X_RIGHT = PAD_HORIZ + half_w + PAD_HORIZ + LABEL_LEFT_OFFSET  # for D

def draw_label(text, x, y):
    draw.text((x, y), text, fill='black', font=label_font)

y = PAD_VERT
# Panel A
draw_label('A   Model construction and evaluation overview', LABEL_X_LEFT, y)
y += LABEL_HEIGHT
canvas.paste(panel_a, (PAD_HORIZ, y))
y += panel_a.height + PAD_VERT

# Panel B
draw_label('B   Predicted-probability distributions on CAMP-Balanced',
           LABEL_X_LEFT, y)
y += LABEL_HEIGHT
canvas.paste(panel_b, (PAD_HORIZ, y))
y += panel_b.height + PAD_VERT

# Panel C + D
draw_label('C   Pairwise AUC comparison', LABEL_X_LEFT, y)
draw_label('D   Top-quartile case enrichment', LABEL_X_RIGHT, y)
y += LABEL_HEIGHT
canvas.paste(panel_c, (PAD_HORIZ, y))
canvas.paste(panel_d, (PAD_HORIZ + half_w + PAD_HORIZ, y))

print(f"\nCombined canvas: {canvas.size}")

# Save as PDF at 300 dpi
canvas.save(OUT_PDF, 'PDF', resolution=300.0)
print(f"Saved: {OUT_PDF}")

# Also save a preview PNG for Quick Look
preview_png = FIGDIR / 'Figure4_preview.png'
canvas.save(preview_png, 'PNG', dpi=(300, 300))
print(f"Saved preview: {preview_png}")
