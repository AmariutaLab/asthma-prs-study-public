"""Generate Figure 4 panels — individual PNG + PDF files for separate review.

This iteration:
  * Panel B — tighter whitespace (smaller hspace, tighter top/bottom margins)
  * Panel C — vertical size halved; in-panel legend text box removed (model
              abbreviations + bootstrap explanation moved to the caption); keep
              colorbar with axis annotation only
  * Panel D — vertical size halved; in-panel legend removed (errorbar / OR=1 /
              bar-color legend moved to caption); per-bar 'OR = X.XX' text removed
              (also covered by caption). Keeps x-axis label + model tick labels.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path
from sklearn.metrics import roc_auc_score

# --- Matplotlib settings for EDITABLE SVG output ---------------------------
# svg.fonttype='none' keeps text as <text> elements (NOT outlined paths), so
# every label remains editable after PowerPoint's "Convert to Shape".
# svg.image_inline=False would write rasterized blobs to external files; we
# keep image_inline=True (the default) but disable rasterization explicitly.
matplotlib.rcParams['svg.fonttype']      = 'none'
matplotlib.rcParams['svg.image_inline']  = True
matplotlib.rcParams['pdf.fonttype']      = 42  # keep PDFs text-editable too
matplotlib.rcParams['ps.fonttype']       = 42

# --- Standardized font sizes for B/C/D (Cell Genomics target ~7-8 pt print
# scale when 16-inch PPT is downsized to 6.85 inches; 12-14 pt at PPT scale).
# Override per-element font sizes via rcParams so individual sns/plt calls
# inherit a consistent look. Per-call fontsize=X args still override these.
matplotlib.rcParams['axes.titlesize']        = 14
matplotlib.rcParams['axes.labelsize']        = 13
matplotlib.rcParams['xtick.labelsize']       = 12
matplotlib.rcParams['ytick.labelsize']       = 12
matplotlib.rcParams['legend.fontsize']       = 11
matplotlib.rcParams['legend.title_fontsize'] = 12

sns.set_theme(style='whitegrid', context='paper', font_scale=1.0)
CB = sns.color_palette('colorblind')

# --- Panel-A-consistent model-class palette (bars only; fonts stay black) ---
CLASS_COLORS = {'PRS':'#7fbf7b','TWAS':'#f0a860','FOCUS':'#e6c84e','UNI':'#a99bd0','XM':'#7aa9db'}
def class_color(key):
    for pre,cls in (('TWAS','TWAS'),('FOCUS','FOCUS'),('PRS_','PRS'),('UNI_','UNI'),('XM_','XM')):
        if key.startswith(pre): return CLASS_COLORS[cls]
    return CLASS_COLORS['TWAS']

ROOT = Path(os.environ.get('FIG4_ROOT', '/Users/nancyh/Desktop/hartwell/gene_model/score/combine'))  # override with FIG4_ROOT env var
PRED = ROOT / 'data' / 'predictions'
FIGDIR = ROOT / 'figures'
FIGDIR.mkdir(parents=True, exist_ok=True)

N_BOOT = 100
BAL_N = 64

# Absolute (in inches) placement of the panel letter + descriptive title so they
# appear at the same x-offset from the left edge of every panel image regardless
# of figure width.
PANEL_LETTER_X_IN  = 0.10
PANEL_TITLE_X_IN   = 0.45
PANEL_TOP_Y_IN     = 0.18
PANEL_LETTER_FS    = 18
PANEL_TITLE_FS     = 15


def panel_label(fig, letter, title):
    fw = fig.get_figwidth()
    fh = fig.get_figheight()
    y_letter = 1.0 - PANEL_TOP_Y_IN / fh
    y_title  = 1.0 - (PANEL_TOP_Y_IN + 0.02) / fh
    fig.text(PANEL_LETTER_X_IN / fw, y_letter, letter,
             fontsize=PANEL_LETTER_FS, fontweight='bold', va='top', ha='left')
    fig.text(PANEL_TITLE_X_IN / fw, y_title, title,
             fontsize=PANEL_TITLE_FS, fontweight='bold', va='top', ha='left')


# ---------------------------------------------------------------- data loaders

def load_consistent(path):
    d = pd.read_csv(path)
    if {'y_pred', 'y_true'} <= set(d.columns):
        return d['y_pred'].astype(float).values, d['y_true'].astype(int).values
    if {'score', 'y_true', 'cohort'} <= set(d.columns):
        d = d[d['cohort'] == 'CAMP_only']
        return d['score'].astype(float).values, d['y_true'].astype(int).values
    raise ValueError(f'unknown format: {path}')


# ---- Dynamic anchor resolvers ----
# The figure previously hard-coded a per-sample prediction path for each
# TWAS-P+T and MA-FOCUS anchor (feature + classifier [+ p-value]). When the
# upstream best-classifier pick changed (e.g. cd4_naive: RF-GridSearch ->
# Gradient Boosting), the hard-coded paths went stale, so Figure 4 and
# Supplementary Table S7 disagreed for the same anchor.
#
# These resolvers pick the classifier (and, for TWAS P+T, the p-value
# threshold) from `best_valid_per_feature.csv` at run time, so Figure 4
# tracks whatever S7 tracks.
import re as _re


def _safe_model_name(m):
    return _re.sub(r'[^A-Za-z0-9_]+', '_', str(m)).strip('_')


def _resolve_predictions(pred_dir, feature, model):
    """Locate the CAMP-only per-sample prediction CSV for
    (feature, model) under `pred_dir/predictions/`.

    Handles safe-name mangling variants ("Ridge (C=0.01)" -> "Ridge_C_0_01" etc.).
    Also falls back to the multi-cohort long-format file (filtering
    cohort == 'CAMP_only' downstream in load_consistent)."""
    safe = _safe_model_name(model)
    candidates = [
        pred_dir / 'predictions' / f'best_consistent__{feature}__{safe}__CAMP_only.csv',
        pred_dir / 'predictions' / f'best_consistent__{feature}__{safe}.csv',
    ]
    for c in candidates:
        if c.exists(): return c
    # Last-ditch glob (handles minor safe-name variants)
    hits = list((pred_dir / 'predictions').glob(f'best_consistent__{feature}__*__CAMP_only.csv'))
    if hits:
        return hits[0]
    raise FileNotFoundError(f'No CAMP-only prediction file for {feature} / {model} under {pred_dir}/predictions/')


def resolve_focus_anchor(mv, feature):
    """Look up the CURRENT best classifier for a MA-FOCUS anchor feature via
    best_valid_per_feature.csv, return the CAMP-only per-sample path."""
    bv = pd.read_csv(PRED / f'meta_model_{mv}' / 'best_valid_per_feature.csv')
    row = bv[bv['Feature'] == feature]
    if row.empty:
        raise ValueError(f'{feature} not in meta_model_{mv}/best_valid_per_feature.csv')
    row = row.iloc[0]
    model = row['Model']
    path = _resolve_predictions(PRED / f'meta_model_{mv}', feature, model)
    print(f'  FOCUS {mv} anchor={feature} -> classifier={model} (AUC={row["CAMP_ONLY_BAL_AUC"]:.4f}) [{path.name}]')
    return path


def resolve_twas_pt_anchor(mv, feature, pvs=('5e-05', '5e-04', '0_005', '0_05')):
    """Sweep all TWAS P+T p-value dirs; pick the pval whose
    best_valid_per_feature.csv gives the HIGHEST CAMP_ONLY_BAL_AUC for
    `feature`. Return the CAMP-only per-sample path."""
    best = None
    for pv in pvs:
        d = PRED / f'meta_model_{mv}__pval-{pv}'
        bv_path = d / 'best_valid_per_feature.csv'
        if not bv_path.exists(): continue
        bv = pd.read_csv(bv_path)
        m = bv[bv['Feature'] == feature]
        if m.empty: continue
        r = m.iloc[0]
        auc = float(r['CAMP_ONLY_BAL_AUC'])
        if best is None or auc > best['auc']:
            best = {'auc': auc, 'pv': pv, 'model': r['Model'], 'dir': d}
    if best is None:
        raise ValueError(f'{feature} not in any meta_model_{mv}__pval-*/best_valid_per_feature.csv')
    path = _resolve_predictions(best['dir'], feature, best['model'])
    print(f'  TWAS P+T {mv} anchor={feature} -> pval={best["pv"]}, '
          f'classifier={best["model"]} (AUC={best["auc"]:.4f}) [{path.name}]')
    return path


def load_prs(method, eval_set='CAMP-only', config='PRS + Ancestry PCs'):
    prs = pd.read_csv(PRED / 'prscs_evaluation' / 'prs_predictions.csv')
    sub = prs[(prs['method'] == method) & (prs['eval_set'] == eval_set) & (prs['config'] == config)]
    return sub['score'].values, sub['y_true'].astype(int).values


def load_unified(mv):
    bal = pd.read_csv(PRED / f'meta_model_{mv}' / 'unified_balanced_results.csv')
    best_method = bal.loc[bal['AUC'].idxmax(), 'Method']
    u = pd.read_csv(PRED / f'meta_model_{mv}' / 'predictions' / 'unified_ptrs_long.csv')
    s = u[(u['method'] == best_method) & (u['cohort'] == 'CAMP_only')]
    return s['score'].values, s['y_true'].astype(int).values


def load_crossmodal(prs_type):
    a = pd.read_csv(PRED / 'integrated_ptrs_prs_combined' / 'all_results.csv')
    sub = a[(a['Approach'] == 'Direct') & (a['Eval_Set'] == 'CAMP-only Balanced') &
            (a['Method'].str.contains(f'\\({prs_type}\\)', regex=True))]
    method_full = sub.loc[sub['AUC'].idxmax(), 'Method']
    classifier = method_full.replace(f'Direct ({prs_type}) + ', '')
    dp = pd.read_csv(PRED / 'integrated_ptrs_prs_combined' / 'direct_predictions.csv')
    s = dp[(dp['method'] == classifier) & (dp['prs_type'] == prs_type) & (dp['eval_set'] == 'CAMP-only')]
    return s['score'].values, s['y_true'].astype(int).values


ut_s, ut_y = load_unified('tissue')
uc_s, uc_y = load_unified('ct')
xm_cs_s, xm_cs_y   = load_crossmodal('PRS-CS')
xm_csx_s, xm_csx_y = load_crossmodal('PRS-CSx')

print("\n---- Resolving anchor classifiers from current best_valid_per_feature.csv ----")
models = [
    ('TWAS Tissue\n(Spleen)',
     *load_consistent(resolve_twas_pt_anchor('tissue', 'Spleen')),
     'r1', 0, 'TWAS_T'),
    ('TWAS Celltype\n' + r'(CD56$^{\mathrm{bright}}$ NK)',
     *load_consistent(resolve_twas_pt_anchor('ct', 'nk_cd56bright')),
     'r1', 1, 'TWAS_CT'),
    ('MA-FOCUS Tissue\n(Esophagus Mucosa)',
     *load_consistent(resolve_focus_anchor('tissue', 'Esophagus_Mucosa')),
     'r1', 2, 'FOCUS_T'),
    ('MA-FOCUS Celltype\n' + r'(CD4$^+$ Naive T)',
     *load_consistent(resolve_focus_anchor('ct', 'cd4_naive')),
     'r1', 3, 'FOCUS_CT'),
    ('PRS-CS\n(default)',     *load_prs('PRS-CS'),     'r2', 0, 'PRS_CS'),
    ('PRS-CSx\n(default)',    *load_prs('PRS-CSx'),    'r2', 1, 'PRS_CSx'),
    ('Unified\nTissue PTRS',  ut_s, ut_y,              'r2', 2, 'UNI_T'),
    ('Unified\nCelltype PTRS', uc_s, uc_y,             'r2', 3, 'UNI_CT'),
    ('Cross-modal\n+ PRS-CS', xm_cs_s, xm_cs_y,        'r2', 4, 'XM_CS'),
    ('Cross-modal\n+ PRS-CSx', xm_csx_s, xm_csx_y,     'r2', 5, 'XM_CSx'),
]
print(f"Loaded {len(models)} models")


def paired_bootstrap(models_list, n_boot=N_BOOT, bal_n=BAL_N):
    n_models = len(models_list)
    aucs = np.full((n_boot, n_models), np.nan)
    ors  = np.full((n_boot, n_models), np.nan)
    for j, m in enumerate(models_list):
        _, scores, y_true, _, _, _ = m
        scores = np.asarray(scores); y_true = np.asarray(y_true).astype(int)
        case_idx = np.where(y_true == 1)[0]
        ctrl_idx = np.where(y_true == 0)[0]
        n_each = min(len(case_idx), len(ctrl_idx), bal_n)
        for i in range(n_boot):
            rng = np.random.RandomState(i)
            cs = rng.choice(case_idx, size=n_each, replace=False)
            ct = rng.choice(ctrl_idx, size=n_each, replace=False) if len(ctrl_idx) > n_each else ctrl_idx
            idx = np.concatenate([cs, ct])
            y_b = y_true[idx]; s_b = scores[idx]
            aucs[i, j] = roc_auc_score(y_b, s_b)
            q1, q4 = np.quantile(s_b, [0.25, 0.75])
            top = y_b[s_b >= q4]; bot = y_b[s_b <= q1]
            a, b = top.sum(), len(top) - top.sum()
            c, d = bot.sum(), len(bot) - bot.sum()
            ors[i, j] = (a * d) / (b * c) if (b > 0 and c > 0) else np.nan
    return aucs, ors


print(f"\nRunning paired bootstrap × {N_BOOT} ...")
aucs, ors = paired_bootstrap(models)
auc_mean = np.nanmean(aucs, axis=0)
auc_std  = np.nanstd(aucs, axis=0)
or_mean  = np.nanmean(ors, axis=0)
or_lo    = np.nanpercentile(ors, 2.5, axis=0)
or_hi    = np.nanpercentile(ors, 97.5, axis=0)

n_models = len(models)
delta_auc = np.zeros((n_models, n_models))
pval = np.ones((n_models, n_models))
for i in range(n_models):
    for j in range(n_models):
        if i == j: continue
        diff = aucs[:, i] - aucs[:, j]
        delta_auc[i, j] = diff.mean()
        pval[i, j] = (diff <= 0).mean()




# ============================================================================ B
from scipy.stats import ttest_ind
from matplotlib.lines import Line2D


def _sig_parts(p):
    # Capital P (Cell Press convention). Italic P will be added at typesetting.
    if np.isnan(p):    return 'n/a', ''
    if p >= 0.05:      return 'ns', ''
    if p < 1e-3:
        s = f'{p:.1e}'.replace('e-0', 'e-').replace('e+0', 'e+')
        return f'P = {s}', '***'
    if p < 0.01:
        return f'P = {p:.3f}', '**'
    return f'P = {p:.3f}', '*'


def violin_row(ax, row_models, fig, show_legend=True, violin_width=0.97,
               standardize='zscore'):
    """Split-violin case-vs-control panel.

    `standardize`:
      'zscore' — Per-model z-standardization: for each model, subtract the
                 pooled (cases + controls) mean and divide by the pooled SD
                 before plotting. Every violin ends up centered at 0 with unit
                 spread, so classifier-specific score compression / expansion
                 no longer distorts the visual comparison. The Welch t-test
                 P-value is scale-invariant, so significance annotations are
                 unchanged; AUC annotations are also unaffected (they come
                 from the paired bootstrap on the raw scores).
      'raw'    — Plot the raw predicted probability of asthma on a shared
                 [0, 1] axis (legacy behavior).
    """
    rows, labels = [], []
    # Visual clip for z-scored violins: a small number of rare outliers
    # (e.g. PRS-CSx has a single sample at z ≈ 6.7) can stretch the KDE tail
    # past the P-value annotation. Winsorize at ±Z_CLIP purely for display —
    # the raw scores are still used for AUC and Welch-P computation below.
    Z_CLIP = 3.5
    for m in row_models:
        label, scores, y_true, _, _, _ = m
        sc = np.asarray(scores, dtype=float)
        if standardize == 'zscore':
            mu = sc.mean()
            sd = sc.std(ddof=0) or 1.0  # guard against constant-prediction models
            sc = (sc - mu) / sd
            sc = np.clip(sc, -Z_CLIP, Z_CLIP)   # visualization-only winsorization
        for s, y in zip(sc, y_true):
            rows.append({'value': float(s), 'class': 'Case' if y == 1 else 'Control', 'model': label})
        labels.append(label)
    df = pd.DataFrame(rows)
    # Palette anchored to Panel A's schematic boxes:
    #   Control (blue) = multimodal (XM) box blue     -> '#7aa9db'
    #   Case (orange)  = single-feature (TWAS/FOCUS)  -> '#fce5cd'  (tan/peach,
    #     sampled from the "Single tissue/cell-type gene-level" box background
    #     in figure4_editable_fix.pptx.png).
    _control_color = CLASS_COLORS['XM']  # '#7aa9db'
    _case_color    = '#fce5cd'           # Panel A single-feature box background
    sns.violinplot(data=df, x='model', y='value', hue='class', split=True, inner='quartile',
                   cut=0, palette={'Control': _control_color, 'Case': _case_color},
                   order=labels, hue_order=['Control', 'Case'], ax=ax,
                   linewidth=0.9, width=violin_width)
    # Always remove the default seaborn legend — we either build our own (split,
    # below the plot) or skip the legend entirely.
    _default_leg = ax.get_legend()
    if _default_leg is not None:
        _default_leg.remove()
    ax.set_xlabel('')
    if standardize == 'zscore':
        ax.set_ylabel('Standardized predicted probability (z)', fontsize=13)
    else:
        ax.set_ylabel('Predicted probability of asthma', fontsize=13)
    ax.set_xlim(-0.55, len(row_models) - 0.45)
    if show_legend:
        # Option B — one horizontal strip with ONLY the visual keys that appear
        # on the plot. Bracket anchor symbol removed: split-violins already
        # visually pair Case (right) vs Control (left), so the bracket is redundant.
        # "Inner lines = Q1/Q2/Q3" and "AUC = mean ± SD across 100 iters" live
        # in the caption.
        strip_handles = [
            mpatches.Patch(facecolor=_control_color, edgecolor='black', linewidth=0.5,
                           label='Control'),
            mpatches.Patch(facecolor=_case_color, edgecolor='black', linewidth=0.5,
                           label='Case'),
            Line2D([], [], color='none', marker='', linestyle='',
                   label='***  p < 0.001'),
            Line2D([], [], color='none', marker='', linestyle='',
                   label='**  p < 0.01'),
            Line2D([], [], color='none', marker='', linestyle='',
                   label='*  p < 0.05'),
            Line2D([], [], color='none', marker='', linestyle='',
                   label='ns  not significant'),
        ]
        fig.legend(
            handles=strip_handles,
            loc='lower center',
            # Lift legend slightly off the very-bottom edge so it doesn't
            # crowd the two-line x-tick labels above it.
            bbox_to_anchor=(0.5, 0.02),
            bbox_transform=fig.transFigure,
            ncol=len(strip_handles),
            fontsize=11, frameon=True, framealpha=0.95,
            handlelength=1.2, handletextpad=0.4, borderpad=0.4,
            columnspacing=1.4,
        )
    else:
        leg = ax.get_legend()
        if leg is not None: leg.remove()

    # Per-panel y-limits + annotation heights. In z-space, ~99% of data
    # sits inside [-3, +3], so we anchor the plot from -3.5 to +3.5 and
    # place AUC/P annotations in a headroom band above the violins.
    if standardize == 'zscore':
        # Two annotation rows share the headroom above the violins:
        #   * P-value / ns row at y_pval — placed with ≥1 z-unit clearance
        #     above the tallest observed violin tail (~z=3.5), so no
        #     collision with the top of the violin.
        #   * AUC bounding-box row at y_auc — placed with ~1 z-unit clearance
        #     above y_pval and ~0.4 z-units below the top of the axis, so no
        #     collision with the P-value row and no clipping at the top.
        ax.set_ylim(-3.6, 7.4)
        ax.set_yticks([-3, -2, -1, 0, 1, 2, 3])
        ax.set_yticklabels(['-3', '-2', '-1', '0', '1', '2', '3'])
        y_auc  = 6.9
        y_pval = 4.7
    else:
        # Raw predicted-probability axis (legacy).
        ax.set_ylim(0.0, 1.18)
        ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1.0'])
        y_auc  = 1.14
        y_pval = 0.95
    for k, m in enumerate(row_models):
        _, scores, y_true, _, _, _ = m
        scores = np.asarray(scores); y_true = np.asarray(y_true).astype(int)
        case_v = scores[y_true == 1]; ctrl_v = scores[y_true == 0]
        if len(case_v) > 0 and len(ctrl_v) > 0:
            # Welch P is scale-invariant — computing it on raw scores is
            # equivalent to computing it on z-scored scores.
            _, p = ttest_ind(case_v, ctrl_v, equal_var=False, nan_policy='omit')
        else:
            p = np.nan
        # Inline p-value + significance marker (no bracket, no separate * line)
        pval_text, star_text = _sig_parts(p)
        inline = f'{pval_text}  {star_text}' if star_text else pval_text
        ax.text(k, y_pval, inline,
                ha='center', va='bottom', fontsize=12, fontweight='normal')
        idx = next(i for i, mm in enumerate(models) if mm[5] == m[5])
        a = auc_mean[idx]; s = auc_std[idx]
        ax.text(k, y_auc, f'AUC {a:.3f} ± {s:.3f}',
                ha='center', va='top', fontsize=12,
                bbox=dict(boxstyle='round,pad=0.10', facecolor='white',
                          edgecolor='lightgray', alpha=0.92))
    # Two-line labels (prefix on top, parenthetical on bottom) — no rotation needed
    plt.setp(ax.get_xticklabels(), rotation=0, ha='center')
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)


# Single-row layout: all 10 violins in one subplot. Wider figure than the 2-row
# layout but much shorter vertically. violin_width=0.80 packs violins close to
# each other while leaving a hairline gap. The two-block legend goes BELOW the
# axes (not to the right), so we use almost the full figure width for the plot.
# Reduce the vertical canvas by ~25% vs the original — the z-standardized
# violins live in a compact [-3, +3] envelope, so a shorter panel works. The
# bottom band is generous enough to fit the two-line x-tick labels AND the
# horizontal legend strip with a visible gap between them.
figB = plt.figure(figsize=(18, 3.7))
# [left, bottom, width, height] — bottom band (~36% of figure) holds:
#   - the axes' x-tick labels (two lines each, ~0.55 in tall)
#   - a clear gap
#   - the horizontal legend strip (~0.30 in tall)
# so 3.7*0.36 = ~1.33 in vertical budget below the axes, comfortably absorbing both.
ax_b = figB.add_axes([0.055, 0.36, 0.93, 0.53])
# Panel B column order (per PI request):
#   PRS-CS, PRS-CSx -> TWAS P+T (T, CT) -> MA-FOCUS (T, CT) -> Unified (T, CT) -> Cross-modal (CS, CSx)
_PANEL_B_ORDER = ['PRS_CS', 'PRS_CSx',
                  'TWAS_T', 'TWAS_CT',
                  'FOCUS_T', 'FOCUS_CT',
                  'UNI_T',   'UNI_CT',
                  'XM_CS',   'XM_CSx']
_model_by_key = {m[5]: m for m in models}
all_violin_models = [_model_by_key[k] for k in _PANEL_B_ORDER]
violin_row(ax_b, all_violin_models, figB, show_legend=True,
           violin_width=0.95, standardize='zscore')
# In-figure title removed for Panel B — the slide-level title in PPT carries it.
# Keep panel_label off for B/C/D so the embedded SVG has only data.
# panel_label(figB, 'B', 'Predicted-probability distributions on CAMP-Balanced')
figB.savefig(FIGDIR / 'figure4_panelB_violins.png', dpi=200, bbox_inches='tight')
figB.savefig(FIGDIR / 'figure4_panelB_violins.pdf', bbox_inches='tight')
plt.close(figB)
print("Saved Panel B")


# ============================================================================ C
# Vertical size cut roughly in half (7.5 -> 4.0). No in-panel legend text box;
# explanation lives in the caption now. Keep colorbar with up-arrow annotation.
figC, axC = plt.subplots(figsize=(7, 4.0))
all_labels = [m[5] for m in models]

# Re-order rows/columns by ASCENDING mean AUC
order_c = np.argsort(auc_mean)
labels_c    = [all_labels[i] for i in order_c]
delta_auc_c = delta_auc[np.ix_(order_c, order_c)]
pval_c      = pval[np.ix_(order_c, order_c)]

mat = np.full((n_models, n_models), np.nan)
for i in range(n_models):
    for j in range(i):
        mat[i, j] = delta_auc_c[i, j]

# square=False so cells stretch horizontally and the heatmap stays short.
# rasterized=False forces each cell to be a separate vector <rect> in the SVG,
# so PowerPoint's "Convert to Shape" decomposes the heatmap into 55 editable
# colored rectangles instead of a single rasterized PNG blob.
hm = sns.heatmap(mat, cmap='Reds', vmin=0.0, vmax=0.15,
                 xticklabels=labels_c, yticklabels=labels_c,
                 cbar_kws={'label': r'$\Delta$AUC', 'shrink': 0.85, 'pad': 0.02},
                 mask=np.triu(np.ones_like(mat, dtype=bool)),
                 ax=axC, linewidths=0.4, linecolor='white', square=False)
# Defensively un-rasterize every QuadMesh/Collection that seaborn created
for coll in axC.collections:
    coll.set_rasterized(False)
for i in range(n_models):
    for j in range(i):
        p = pval_c[i, j]
        star = '**' if p < 0.01 else ('*' if p < 0.05 else '')
        if star:
            axC.text(j + 0.5, i + 0.5, star, ha='center', va='center',
                     fontsize=12, color='black', fontweight='bold')
axC.set_title('')
axC.set_xticklabels(labels_c, rotation=45, ha='right', fontsize=12)
axC.set_yticklabels(labels_c, rotation=0, fontsize=12)
axC.set_xlabel('')   # no axis title — captions describe the convention
axC.set_ylabel('')
axC.grid(False)
cbar = axC.collections[0].colorbar
cbar.ax.tick_params(labelsize=9)
cbar.set_label(r'$\Delta$AUC', fontsize=13)
cbar.ax.tick_params(labelsize=11)

# In-panel legend: ONLY the model-class abbreviation key. Everything else
# (bootstrap convention, sort note, * / ** thresholds) lives in the caption.
axC.text(0.98, 0.95,
         "Model class:\n"
         "  XM    = Cross-modal\n"
         "  FOCUS = MA-FOCUS PTRS\n"
         "  UNI   = Unified PTRS\n"
         "  TWAS  = TWAS P+T\n"
         "  PRS   = PRS only",
         transform=axC.transAxes,
         ha='right', va='top', fontsize=11,
         family='monospace',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                   edgecolor='lightgray', alpha=0.95))

figC.tight_layout(rect=[0, 0, 1, 0.88])
# panel_label(figC, 'C', 'Pairwise AUC comparison')   # slide-level title only
figC.savefig(FIGDIR / 'figure4_panelC_pairwise.png', dpi=200, bbox_inches='tight')
figC.savefig(FIGDIR / 'figure4_panelC_pairwise.pdf', bbox_inches='tight')
plt.close(figC)
print("Saved Panel C")


# ============================================================================ D
# Vertical size halved (7.5 -> 4.0). In-panel legend removed (bar-class colors,
# errorbar handle, OR=1 line) — all moved to the caption. Per-bar 'OR = X.XX'
# annotations also removed (caption now lists the headline OR values).
figD, axD = plt.subplots(figsize=(7, 4.0))

order = np.argsort(-or_mean)
all_labels_sorted = [all_labels[i] for i in order]
or_mean_sorted = or_mean[order]
or_lo_sorted   = or_lo[order]
or_hi_sorted   = or_hi[order]

class_colors_sorted = []
for i in order:
    key = models[i][5]
    class_colors_sorted.append(class_color(key))

y_pos = np.arange(n_models)
or_his_clipped = np.minimum(or_hi_sorted, 15)
# Thinner bars so the panel reads as horizontal lines, not big rectangles
BAR_HEIGHT = 0.55
axD.barh(y_pos, or_mean_sorted, height=BAR_HEIGHT,
         xerr=[or_mean_sorted - or_lo_sorted, or_his_clipped - or_mean_sorted],
         color=class_colors_sorted, alpha=0.85, edgecolor='black', linewidth=0.5,
         error_kw={'elinewidth': 1.0, 'capsize': 2.5, 'ecolor': 'black'})
axD.axvline(1.0, color='gray', linestyle='--', linewidth=1.0, alpha=0.7)
axD.set_yticks(y_pos)
axD.set_yticklabels(all_labels_sorted, fontsize=12)
axD.set_xlabel('Odds ratio (top quartile vs bottom quartile of predicted score)', fontsize=13)
axD.set_title('')
axD.invert_yaxis()
axD.grid(axis='x', alpha=0.3)
axD.tick_params(axis='x', labelsize=11)
axD.set_xlim(0, 8)

figD.tight_layout(rect=[0, 0, 1, 0.88])
# panel_label(figD, 'D', 'Top-quartile case enrichment')   # slide-level title only
figD.savefig(FIGDIR / 'figure4_panelD_or.png', dpi=200, bbox_inches='tight')
figD.savefig(FIGDIR / 'figure4_panelD_or.pdf', bbox_inches='tight')
plt.close(figD)
print("Saved Panel D")


# ---- Backing data ----
out = pd.DataFrame({
    'model_key': [m[5] for m in models],
    'model_label': [m[0].replace('\n', ' ') for m in models],
    'row': [m[3] for m in models],
    'pos': [m[4] for m in models],
    'AUC_mean': auc_mean, 'AUC_std': auc_std,
    'OR_mean': or_mean, 'OR_lo95': or_lo, 'OR_hi95': or_hi,
})
out.to_csv(FIGDIR / 'figure4_revised_data.csv', index=False)
print(f"\nAll 4 panels saved in {FIGDIR}")
