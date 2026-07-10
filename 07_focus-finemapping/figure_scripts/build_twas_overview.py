"""
TWAS overview heatmap (supporting panel for the methods / Figure 2 area).

Shows the top 5 canonical asthma-effector genes per modality
(5 for GTEx tissues, 5 for OneK1K cell types) with their TWAS Z-scores
across the FOCUS credible sets they appear in.

Visual encoding:
    Cell color    signed FUSION marginal TWAS Z (RdBu_r, clipped at ±12) —
                  the .dat TWAS.Z, so colour and the * FDR overlay are one
                  self-consistent statistic (NOT the FOCUS fine-mapping z)
    *  marker     FDR < 0.05 (Benjamini-Hochberg across all gene-context
                  TWAS tests within the modality — the FULL test universe,
                  not just credible-set entries)
    ●  marker     FOCUS PIP ≥ 0.8 (high-confidence fine-mapping)
    light gray    gene absent from this context's credible set

Outputs (see CONFIG block):
    combine/figures/twas_overview_heatmap.png  (200 dpi)
    combine/figures/twas_overview_heatmap.pdf  (vector, editable text)
"""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import seaborn as sns
from scipy.stats import false_discovery_control   # requires scipy >= 1.11

# Editable text in SVG / PDF — keep glyph data, don't outline to paths
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42

# =============================================================================
# CONFIG — edit these paths for your environment (local laptop vs server)
# =============================================================================
# Defaults assume the local laptop layout. On the server, override either by
# editing this block or by setting the corresponding environment variables.
ROOT = Path(os.environ.get('HARTWELL_ROOT',
                           '/Users/nancyh/Desktop/hartwell'))
FOCUS_DIR = Path(os.environ.get('FOCUS_DIR', ROOT / 'FOCUS'))

# FOCUS credible-set tables (one row per gene-context pair in the cred set).
# Tissue: prefer the current (non-PIP-thresholded) table; fall back to the
# older PIP-named file only if the current one is absent.
_T_CANDIDATES = ['focus_credset_gene_tissue_39.tsv',
                 'focus_credset_pip_gene_tissue_39.tsv']
T_FILE   = next((FOCUS_DIR / n for n in _T_CANDIDATES if (FOCUS_DIR / n).exists()),
                FOCUS_DIR / _T_CANDIDATES[0])
C_FILE   = FOCUS_DIR / 'focus_credset_gene_tissue_17CT.tsv'
ANN_FILE = FOCUS_DIR / 'gene_annotation.txt.gz'

# Directory containing FUSION TWAS .dat output files, named
# `twas_<context>_chr<N>.dat`. Used for the FDR overlay; if not present,
# the heatmap renders without `*` FDR markers and prints a warning.
# Override via env var TWAS_DAT_DIR on the server.
TWAS_DAT_DIR_T = Path(os.environ.get('TWAS_DAT_DIR_T',
                                     ROOT / 'TWAS' / 'TWAS_results'))
TWAS_DAT_DIR_C = Path(os.environ.get('TWAS_DAT_DIR_C',
                                     ROOT / 'TWAS' / 'TWAS_results_ct'))

# Output directory for the figure
OUTDIR = Path(os.environ.get('FIG_OUTDIR',
                             ROOT / 'gene_model' / 'score' / 'combine' / 'figures'))
OUTDIR.mkdir(parents=True, exist_ok=True)

# Cache for the (slow) .dat aggregation. The per-modality aggregated table
# [gene, context, p, z] is pickled here and reused on subsequent runs whenever
# the underlying .dat files are unchanged (same set of files, mtimes and sizes).
# Force a rebuild with `--refresh` or by setting TWAS_REFRESH=1.
CACHE_DIR = Path(os.environ.get('TWAS_CACHE_DIR', OUTDIR / '.twas_cache'))

import argparse
_ap = argparse.ArgumentParser(add_help=True)
_ap.add_argument('--refresh', action='store_true',
                 help='Ignore the .dat aggregation cache and re-read every file.')
_args, _ = _ap.parse_known_args()
REFRESH = _args.refresh or os.environ.get('TWAS_REFRESH', '') not in ('', '0', 'false', 'False')

FDR_ALPHA = 0.05  # BH-FDR threshold for the `*` overlay


# =============================================================================
# Helpers
# =============================================================================
def benjamini_hochberg(p):
    """Return BH-FDR q-values. Thin wrapper around scipy.stats so this script
    uses the SAME implementation as the manuscript's reference notebook
    (01b_MFR_SUPPLEMENT_twas-tissue-celltype-miami-plot_FDR-GO-Enrichment.ipynb).
    Inputs are expected to be NA-free (callers filter beforehand)."""
    p = np.asarray(p, dtype=float)
    if len(p) == 0:
        return p
    return false_discovery_control(p, method='bh')


def load_focus_credset(path, ens2hgnc):
    """Load a FOCUS credible-set table; return DataFrame with columns
    [gene, tissues, twas_z, pip]."""
    d = pd.read_csv(path, sep='\t')
    d = d[d['in_cred_set_pop1'] == 1].copy()
    d['twas_z'] = pd.to_numeric(d['twas_z_pop1'], errors='coerce')
    d['pip']    = pd.to_numeric(d['pips_pop1'],   errors='coerce')
    d['gene']   = d['symbol'].map(ens2hgnc).fillna(d['symbol'])
    return d.dropna(subset=['twas_z'])


def aggregate_twas_dat(dat_dir):
    """Walk a directory of FUSION `twas_<context>_chr<N>.dat` files and
    aggregate into a single DataFrame with columns [gene, context, p].
    Returns an empty DataFrame (and prints a warning) if dat_dir is missing
    or empty, so the rest of the script keeps working.
    """
    dat_dir = Path(dat_dir)
    if not dat_dir.exists():
        print(f"  [FDR] TWAS dat dir not found: {dat_dir} — skipping FDR overlay")
        return pd.DataFrame(columns=['gene', 'context', 'p', 'z'])
    # Tissue FUSION outputs are named twas_<context>_chr<N>.dat; cell-type
    # outputs are named TWAS_metaGWAS_<context>_chr<N>.dat. Match both and
    # strip whichever prefix is present to recover the bare context name.
    files = sorted(dat_dir.glob('twas_*_chr*.dat')) + \
            sorted(dat_dir.glob('TWAS_metaGWAS_*_chr*.dat'))

    # --- cache lookup -------------------------------------------------------
    # Signature = the (name, mtime_ns, size) of every input file. If it matches
    # the cached signature we reuse the pickled aggregate and skip the re-read.
    import hashlib, json
    sig = [[f.name, f.stat().st_mtime_ns, f.stat().st_size] for f in files]
    key = hashlib.md5(str(dat_dir.resolve()).encode()).hexdigest()[:12]
    cache_pkl  = CACHE_DIR / f'twas_agg_{key}.pkl'
    cache_meta = CACHE_DIR / f'twas_agg_{key}.json'
    if not REFRESH and cache_pkl.exists() and cache_meta.exists():
        try:
            if json.loads(cache_meta.read_text()).get('sig') == sig:
                cached = pd.read_pickle(cache_pkl)
                print(f"  [cache] {len(cached):,} (gene,context) rows reused "
                      f"from {cache_pkl.name} ({len(files)} .dat files)")
                return cached
            print("  [cache] .dat files changed — rebuilding aggregate")
        except Exception as e:
            print(f"  [cache] ignoring unreadable cache ({e}) — rebuilding")

    rows = []
    for f in files:
        ctx = f.stem
        for pref in ('TWAS_metaGWAS_', 'twas_'):
            if ctx.startswith(pref):
                ctx = ctx[len(pref):]
                break
        ctx = ctx.rsplit('_chr', 1)[0]
        try:
            d = pd.read_csv(f, sep=r'\s+', engine='python')
        except Exception as e:
            print(f"  [FDR] WARNING: failed to read {f.name}: {e}")
            continue
        if 'ID' not in d.columns or 'TWAS.P' not in d.columns:
            continue
        cols = {'ID': 'gene', 'TWAS.P': 'p'}
        if 'TWAS.Z' in d.columns:
            cols['TWAS.Z'] = 'z'
        d = d[list(cols)].rename(columns=cols)
        d['context'] = ctx
        d['p'] = pd.to_numeric(d['p'], errors='coerce')
        d['z'] = pd.to_numeric(d['z'], errors='coerce') if 'z' in d.columns else np.nan
        d = d.dropna(subset=['p'])
        rows.append(d[['gene', 'context', 'p', 'z']])
    if not rows:
        print(f"  [FDR] No usable TWAS .dat files in {dat_dir}")
        return pd.DataFrame(columns=['gene', 'context', 'p', 'z'])
    out = pd.concat(rows, ignore_index=True)
    # If a (gene, context) pair appears more than once (rare — same gene tested
    # by multiple sub-models), keep the smallest p so FDR is conservative; the
    # retained row's TWAS.Z is the marginal z used for the heatmap colour.
    out = out.sort_values('p').drop_duplicates(['gene', 'context'])
    # --- cache write --------------------------------------------------------
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        out.to_pickle(cache_pkl)
        cache_meta.write_text(json.dumps({'sig': sig, 'nfiles': len(files)}))
        print(f"  [cache] wrote {len(out):,} rows to {cache_pkl.name}")
    except Exception as e:
        print(f"  [cache] WARNING: could not write cache ({e})")
    return out


def compute_fdr_sig_set(twas_df, alpha=0.05):
    """Apply BH-FDR within modality; return set of (gene, context) tuples
    with q < alpha. Empty set if twas_df is empty."""
    if twas_df.empty:
        return set()
    q = benjamini_hochberg(twas_df['p'].values)
    twas_df = twas_df.assign(q=q)
    sig = twas_df[twas_df['q'] < alpha]
    print(f"  [FDR] {len(sig):,} / {len(twas_df):,} gene-context pairs FDR<{alpha}")
    return set(zip(sig['gene'].astype(str), sig['context'].astype(str)))


# Pretty column labels (mathtext where helpful, shortened tissue names)
def pretty(s):
    # T-cell subset abbreviations kept fully upper-case (not .title()-cased)
    SUBSET = {'ctl': 'CTL', 'tcm': 'TCM', 'tem': 'TEM'}
    if s.startswith('cd4_'):
        sub = s.split('cd4_', 1)[1]
        return r'CD4$^+$ ' + SUBSET.get(sub, sub.replace('_', ' ').title())
    if s.startswith('cd8_'):
        sub = s.split('cd8_', 1)[1]
        return r'CD8$^+$ ' + SUBSET.get(sub, sub.replace('_', ' ').title())
    if s == 'nk':                return 'NK'
    if s == 'nk_cd56bright':     return r'NK CD56$^{bright}$'
    if s == 'treg':              return 'Treg'
    if s == 'mait':              return 'MAIT'
    if s == 'b_naive':           return 'B naive'
    if s == 'b_intermediate':    return 'B intermediate'
    if s == 'b_memory':          return 'B memory'
    if s == 'dc':                return 'DC'
    if s == 'plasmablast':       return 'Plasmablast'
    if s == 'monocyte':          return 'Monocyte'
    if s == 'cd14_monocyte':     return 'CD14 mono'
    SHORT = {
        'adipose_subcutaneous':            'Adipose (subcut.)',
        'adipose_visceral_omentum':        'Adipose (visc.)',
        'artery_aorta':                    'Artery (aorta)',
        'artery_coronary':                 'Artery (coronary)',
        'artery_tibial':                   'Artery (tibial)',
        'brain_substantia_nigra':          'Brain (subst. nigra)',
        'brain_basalganglia':              'Brain (basal gang.)',
        'cells_cultured_fibroblasts':      'Fibroblasts (cult.)',
        'colon_transverse':                'Colon (transv.)',
        'colon_sigmoid':                   'Colon (sigm.)',
        'esophagus_mucosa':                'Esophagus mucosa',
        'esophagus_muscularis':            'Esophagus muscul.',
        'heart_atrial_appendage':          'Heart (atrial app.)',
        'skin_not_sun_exposed_suprapubic': 'Skin (not sun-exp.)',
        'skin_sun_exposed_lower_leg':      'Skin (sun-exp.)',
        'breast_mammary_tissue':           'Breast',
        'whole_blood':                     'Whole blood',
        'stomach':                         'Stomach',
        'testis':                          'Testis',
        'thyroid':                         'Thyroid',
        'spleen':                          'Spleen',
        'pancreas':                        'Pancreas',
        'lung':                            'Lung',
    }
    if s in SHORT:
        return SHORT[s]
    return s.replace('_', ' ').title()


# =============================================================================
# Load Ensembl → HGNC + FOCUS credible sets
# =============================================================================
ann = pd.read_csv(ANN_FILE, sep='\t', header=None, compression='gzip',
                  usecols=[3, 6], names=['hgnc', 'ens_versioned'])
ann['ens'] = ann['ens_versioned'].str.split('.').str[0]
ens2hgnc = dict(zip(ann['ens'], ann['hgnc']))

dt = load_focus_credset(T_FILE, ens2hgnc)   # gene-tissue credible set
dc = load_focus_credset(C_FILE, ens2hgnc)   # gene-celltype credible set


# =============================================================================
# Load + FDR-correct full TWAS test universe (per modality)
# =============================================================================
print("Aggregating TWAS .dat files for FDR overlay...")
print(f"  Tissue dir:    {TWAS_DAT_DIR_T}")
twas_t = aggregate_twas_dat(TWAS_DAT_DIR_T)
fdr_sig_t = compute_fdr_sig_set(twas_t, alpha=FDR_ALPHA)

print(f"  Cell-type dir: {TWAS_DAT_DIR_C}")
twas_c = aggregate_twas_dat(TWAS_DAT_DIR_C)
fdr_sig_c = compute_fdr_sig_set(twas_c, alpha=FDR_ALPHA)


# =============================================================================
# Top-5 asthma effector genes per panel
# Picked for (a) appearing in FOCUS credible set, (b) being a named asthma
# gene a reader will recognize, (c) having strong FUSION TWAS signal.
# Cell-type panel chosen for lineage spread (CD4 / NK / CD8) and to avoid
# genes whose signal was inflated only under the FOCUS z (TAPBP, HLA-DPA1).
# =============================================================================
GENES_T = ['TSLP', 'STAT6', 'IL1RL1', 'GATA3', 'GSDMB']           # tissue panel
GENES_C = ['HLA-DPB1', 'IRF1', 'TNFAIP3', 'IL18R1', 'CD247']      # cell-type panel

print()
print(f"Tissue-panel genes:    {GENES_T}")
print(f"Cell-type-panel genes: {GENES_C}")


def make_matrix(d, gene_order):
    """Pivot gene × context for the shortlist. Drop columns with no entries
    among the shown genes. Returns (Z_pivot, P_pivot, raw_cols) where
    raw_cols is the unpretty-formatted column list for FDR lookup."""
    d2 = d[d['gene'].isin(gene_order)].copy()
    d2['abs_z'] = d2['twas_z'].abs()
    d2 = d2.sort_values('abs_z', ascending=False).drop_duplicates(['gene', 'tissues'])
    pivot_z = d2.pivot(index='gene', columns='tissues', values='twas_z').reindex(gene_order)
    pivot_p = d2.pivot(index='gene', columns='tissues', values='pip').reindex(gene_order)
    keep = pivot_z.columns[pivot_z.notna().any(axis=0)].tolist()
    col_counts = pivot_z[keep].notna().sum()
    keep = sorted(keep, key=lambda c: (-col_counts[c], c))
    return pivot_z[keep], pivot_p[keep], list(keep)


Z_t, P_t, raw_t = make_matrix(dt, GENES_T)
Z_c, P_c, raw_c = make_matrix(dc, GENES_C)


# Cell COLOUR = FUSION marginal TWAS.Z (from the .dat files), not the FOCUS
# fine-mapping twas_z. The credible-set table still decides WHICH cells are
# shown (and supplies the ● PIP markers); we just recolour each shown cell by
# its FUSION TWAS.Z, matched on (gene symbol, context) — the .dat ID column is
# the HGNC symbol in both modalities. This makes colour and the * FDR overlay
# one self-consistent FUSION statistic.
def fusion_z_matrix(Z_focus, raw_cols, twas_df):
    zmap = {(str(g), str(c)): z for g, c, z in
            zip(twas_df['gene'], twas_df['context'], twas_df['z'])}
    Zf = Z_focus.copy()
    for i, gene in enumerate(Zf.index):
        for j, ctx in enumerate(raw_cols):
            if pd.isna(Z_focus.iloc[i, j]):
                continue                       # not in credible set -> stays gray
            Zf.iloc[i, j] = zmap.get((str(gene), str(ctx)), np.nan)
    return Zf


Z_t = fusion_z_matrix(Z_t, raw_t, twas_t)
Z_c = fusion_z_matrix(Z_c, raw_c, twas_c)


# Append a dummy '…' column + '…' row to visually signal "more genes / more
# contexts shown elsewhere". The dummy cells are NaN (so they render as the
# absent-gene gray), and draw() overlays a '…' glyph in each one.
ELLIPSIS = '…'

def add_ellipsis_extension(Z, P, raw_cols):
    Z2 = Z.copy()
    P2 = P.copy()
    raw2 = list(raw_cols)
    Z2[ELLIPSIS] = np.nan          # dummy column on the right
    P2[ELLIPSIS] = np.nan
    raw2.append(ELLIPSIS)
    Z2.loc[ELLIPSIS] = np.nan      # dummy row at the bottom
    P2.loc[ELLIPSIS] = np.nan
    return Z2, P2, raw2


Z_t, P_t, raw_t = add_ellipsis_extension(Z_t, P_t, raw_t)
Z_c, P_c, raw_c = add_ellipsis_extension(Z_c, P_c, raw_c)


# Pretty-rename for display only; keep raw_t / raw_c for FDR lookup
Z_t.columns = [pretty(c) for c in Z_t.columns]
P_t.columns = [pretty(c) for c in P_t.columns]
Z_c.columns = [pretty(c) for c in Z_c.columns]
P_c.columns = [pretty(c) for c in P_c.columns]


# =============================================================================
# Plot
# =============================================================================
sns.set_theme(style='white', context='paper')

n_t = Z_t.shape[1]
n_c = Z_c.shape[1]
CELL_W = 0.62          # horizontal cell size (narrower)
ROW_H  = 0.52          # vertical cell size (shorter -> landscape cells)
w_left  = max(3.0, n_t * CELL_W)
w_right = max(3.0, n_c * CELL_W)
fig_h   = max(4.0, len(GENES_T) * ROW_H + 2.6)
fig = plt.figure(figsize=(w_left + w_right + 2.2, fig_h))
gs = fig.add_gridspec(1, 3,
                      width_ratios=[w_left, w_right, 0.30],
                      wspace=0.375,  # gap between the two panels (halved)
                      left=0.09, right=0.94, top=0.82, bottom=0.46)
ax_t = fig.add_subplot(gs[0, 0])
ax_c = fig.add_subplot(gs[0, 1])
cax  = fig.add_subplot(gs[0, 2])

ZMAX = 12.0

# Lightened RdBu_r: blend every colour toward white so the heatmap and its
# colour bar read softer/paler (LIGHTEN is the fraction pulled toward white).
LIGHTEN = 0.45
_rdbu = plt.cm.RdBu_r(np.linspace(0, 1, 256))
_rdbu[:, :3] = _rdbu[:, :3] * (1 - LIGHTEN) + LIGHTEN
CMAP = mcolors.ListedColormap(_rdbu)


def draw(ax, Z, P, raw_cols, fdr_sig, panel_label=None):
    arr = Z.values.copy()
    sns.heatmap(arr, ax=ax, cmap=CMAP, center=0, vmin=-ZMAX, vmax=ZMAX,
                cbar=False, linewidths=0.6, linecolor='white',
                xticklabels=Z.columns, yticklabels=Z.index,
                mask=np.isnan(arr),
                square=False)   # cells fill the axes box: wide (CELL_W) & short (ROW_H)

    # Index of the appended dummy '…' row / column (always the last one).
    n_rows, n_cols = arr.shape
    ell_row = n_rows - 1   # bottom row is the '…' row
    ell_col = n_cols - 1   # right-most col is the '…' col

    # Layered overlays. When a cell has both, place * (FDR) and ● (PIP) side by
    # side horizontally (left / right of center) so they don't collide.
    # Skip the dummy '…' row/col (handled below).
    for i, gene in enumerate(Z.index):
        if i == ell_row:
            continue
        for j, raw_ctx in enumerate(raw_cols):
            if j == ell_col:
                continue
            p = P.iloc[i, j]
            if pd.isna(p):
                continue
            is_fdr = (gene, raw_ctx) in fdr_sig
            is_pip = p >= 0.8
            if is_fdr and is_pip:
                # both — sit horizontally: * on the left, ● on the right
                ax.text(j + 0.35, i + 0.5, '*', ha='center', va='center',
                        fontsize=15, color='black', fontweight='bold')
                ax.text(j + 0.66, i + 0.5, '●', ha='center', va='center',
                        fontsize=8, color='black')
            elif is_fdr:
                ax.text(j + 0.5, i + 0.5, '*', ha='center', va='center',
                        fontsize=15, color='black', fontweight='bold')
            elif is_pip:
                ax.text(j + 0.5, i + 0.5, '●', ha='center', va='center',
                        fontsize=9, color='black')

    # The dummy '…' row / column cells are left empty (blank); the '…' appears
    # only as the last x-axis and y-axis tick label (from the '…' index/column).
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right',
                       rotation_mode='anchor', fontsize=10)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=11,
                       style='italic', fontweight='bold')
    ax.set_xlabel('')
    if panel_label:
        ax.set_title(panel_label, fontsize=12, fontweight='bold', loc='left', pad=6)
    ax.set_facecolor('#f5f5f5')


draw(ax_t, Z_t, P_t, raw_t, fdr_sig_t, panel_label='Bulk tissues (GTEx)')
draw(ax_c, Z_c, P_c, raw_c, fdr_sig_c, panel_label='PBMC cell types (OneK1K)')

# Shared colorbar
norm = mcolors.Normalize(vmin=-ZMAX, vmax=ZMAX)
cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP),
                    cax=cax, orientation='vertical')
cbar.set_label('TWAS Z-score', fontsize=9)
cbar.ax.tick_params(labelsize=8)

# Figure-level title.
# The FOCUS PIP dot in the legend is drawn as its OWN text piece so we can size
# it (slightly larger than a mathtext $\bullet$) and raise it (centred with the
# caps) independently of the surrounding words. The title is laid out as three
# left-anchored pieces — [left text][dot][right text] — measuring each piece's
# width with the renderer so they butt together seamlessly.
TITLE_Y     = 0.965   # baseline of the title text (figure coords)
DOT_SIZE    = 11      # dot point size (a touch larger than the 12-pt words' $\bullet$)
_left_txt   = ('TWAS associations for top asthma effector genes  '
               '(* FDR < 0.05,  ')
_right_txt  = r' FOCUS PIP $\geq$ 0.8)'

fig.canvas.draw()                       # realise a renderer for width measurement
_rend = fig.canvas.get_renderer()
_inv  = fig.transFigure.inverted()

def _place(x, txt, size, va='baseline', y=TITLE_Y):
    t = fig.text(x, y, txt, fontsize=size, fontweight='bold', ha='left', va=va)
    bb = t.get_window_extent(renderer=_rend)
    (x1, _), (_, ytop) = _inv.transform(bb.max), _inv.transform(bb.max)
    return t, x1

# Left words (baseline at TITLE_Y). Measure their cap top so the dot can be
# centred on the caps rather than on the baseline (a baseline-centred dot reads
# far too low next to capital letters).
_x = 0.10
_tleft, _x = _place(_x, _left_txt, 12)
_bb = _tleft.get_window_extent(renderer=_rend)
_ytop = _inv.transform((_bb.x0, _bb.y1))[1]        # cap/ascender top in fig coords
_dot_y = (TITLE_Y + _ytop) / 2.0                   # midway baseline↔cap-top
_tdot, _x = _place(_x, '●', DOT_SIZE, va='center', y=_dot_y)
_tright, _x = _place(_x, _right_txt, 12)


out_png = OUTDIR / 'twas_overview_heatmap.png'
out_pdf = OUTDIR / 'twas_overview_heatmap.pdf'
fig.savefig(out_png, dpi=200, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')
print()
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
print()
print(f"Tissue-panel genes ({len(GENES_T)}): {GENES_T}")
print(f"  → tissues shown: {raw_t}")
print(f"Cell-type-panel genes ({len(GENES_C)}): {GENES_C}")
print(f"  → cell types shown: {raw_c}")
