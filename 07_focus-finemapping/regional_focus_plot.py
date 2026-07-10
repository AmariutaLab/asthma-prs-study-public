"""
Gene-level regional FOCUS / TWAS plot (LocusZoom-style; each dot is a gene-model).

Two stages, decoupled so the plotted data is frozen and committed:
  1. build_dataset(): raw per-block ma-FOCUS output  ->  one tidy CSV for the panel
  2. regional_plot(): the committed CSV              ->  the figure (PNG + PDF)

Panel C: GTEx bulk tissues, asthma 17q21 locus, focal gene GSDMB.
Panel D (cell types): identical, once raw OneK1K per-block FOCUS output is staged
  (e.g. results_ct/<cell_type>/*.focus.tsv). See --help.

Per dot:  x = gene center (Mb), y = FUSION marginal TWAS Z-score
          (the .dat TWAS.Z, via --fusion-dat-dir), color = FOCUS PIP.
Rings:    red = focal-gene model (one per group);  black = in FOCUS 90% credible set.

Usage
-----
# build the committed dataset, then plot (Panel C, tissues):
python regional_focus_plot.py --focus-dir . --focal-ens ENSG00000073605 \
       --focal-name GSDMB --panel-letter C --group-label Tissues \
       --results-glob 'results/*/*.focus.tsv' --window-mb 1 \
       --csv panel_C_data.csv --out panel_C_17q21_tissues.png

# Panel D on the server (cell types), once raw CT blocks are staged:
python regional_focus_plot.py --focus-dir . --focal-ens <ENSG> \
       --focal-name <GENE> --panel-letter D --group-label 'Cell types' \
       --results-glob 'results_ct/*/*.focus.tsv' --window-mb 1 \
       --csv panel_D_data.csv --out panel_D_celltypes.png

# replot from an already-built CSV (no raw data needed):
python regional_focus_plot.py --csv panel_C_data.csv --plot-only \
       --focal-ens ENSG00000073605 --focal-name GSDMB \
       --panel-letter C --group-label Tissues --out panel_C_17q21_tissues.png
"""
import os
import glob
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# columns kept in the committed CSV
CSV_COLS = ['chrom', 'ensembl', 'symbol', 'group', 'tx_start', 'tx_stop',
            'gene_center', 'twas_z', 'pip', 'in_cred_set', 'x_plot']


# ----------------------------------------------------------------------------
# stage 1: raw ma-FOCUS block output -> committed panel CSV
# ----------------------------------------------------------------------------
def _load_annotation(annotation_path):
    ann = pd.read_csv(annotation_path, sep='\t', header=None, usecols=[3, 6],
                      names=['name', 'ens'], compression='gzip')
    ann['base'] = ann['ens'].astype(str).str.split('.').str[0]
    return dict(zip(ann['base'], ann['name']))


def _fusion_z_map(dat_dir, chrom):
    """Map gene -> FUSION marginal TWAS.Z for one chromosome, keeping the
    smallest-p model per key. Returns (by_ensembl, by_symbol) dicts keyed by
    (id, context). Tissue .dat name their weights by Ensembl id (in the FILE
    column); cell-type .dat name them by symbol (the ID column) — so we index
    both ways and let the caller match on whichever it has."""
    by_ens, by_sym = {}, {}
    files = (sorted(glob.glob(os.path.join(dat_dir, f'twas_*_chr{chrom}.dat'))) +
             sorted(glob.glob(os.path.join(dat_dir, f'TWAS_metaGWAS_*_chr{chrom}.dat'))))
    for f in files:
        base = os.path.basename(f)
        ctx = None
        for pref in ('TWAS_metaGWAS_', 'twas_'):
            if base.startswith(pref):
                ctx = base[len(pref):].rsplit(f'_chr{chrom}.dat', 1)[0]
                break
        if ctx is None:
            continue
        try:
            d = pd.read_csv(f, sep=r'\s+', engine='python')
        except Exception:
            continue
        if 'ID' not in d.columns or 'TWAS.Z' not in d.columns:
            continue
        for _, r in d.iterrows():
            z = pd.to_numeric(r.get('TWAS.Z'), errors='coerce')
            pp = pd.to_numeric(r.get('TWAS.P'), errors='coerce')
            if pd.isna(z):
                continue
            pp = pp if pd.notna(pp) else 1.0
            sym = str(r['ID'])
            mm = re.search(r'(ENSG\d+)', str(r.get('FILE', '')))
            keys = [((sym, ctx), by_sym)]
            if mm:
                keys.append(((mm.group(1), ctx), by_ens))
            for key, dct in keys:
                if key not in dct or pp < dct[key][1]:
                    dct[key] = (z, pp)
    return ({k: v[0] for k, v in by_ens.items()},
            {k: v[0] for k, v in by_sym.items()})


def _find_block(focal_ens, results_glob, focus_dir):
    """Return the raw FOCUS block .tsv that contains the focal gene (all groups)."""
    for f in glob.glob(os.path.join(focus_dir, results_glob)):
        if os.path.getsize(f) == 0:
            continue
        try:
            ids = pd.read_csv(f, sep='\t', usecols=['ens_gene_id'])
        except Exception:
            continue
        if ids['ens_gene_id'].astype(str).str.contains(focal_ens).any():
            return f
    raise FileNotFoundError(f"no block under {results_glob!r} contains {focal_ens}")


def build_dataset(focus_dir, focal_ens, results_glob, window_mb, out_csv, annotation_path,
                  fusion_dat_dir=None):
    """Slice the focal-gene block to a window and write the tidy, frozen CSV."""
    e2n = _load_annotation(annotation_path)
    block = _find_block(focal_ens, results_glob, focus_dir)
    df = pd.read_csv(block, sep='\t')

    df = df[df['ens_gene_id'].astype(str) != 'NULL.MODEL'].copy()
    df['ensembl'] = df['ens_gene_id'].astype(str).str.split('_').str[0]
    df['group'] = df['ens_gene_id'].astype(str).str.split('_', n=1).str[1]   # tissue / cell type
    df['symbol'] = df['ensembl'].map(e2n).fillna(df['ensembl'])
    for c in ['tx_start', 'tx_stop', 'twas_z_pop1', 'pips_pop1']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    df['gene_center'] = (df['tx_start'] + df['tx_stop']) / 2
    df = df.dropna(subset=['gene_center', 'twas_z_pop1'])

    # window around the focal gene's center
    center = df.loc[df['ensembl'] == focal_ens, 'gene_center'].median()
    half = window_mb * 1e6 / 2
    win = df[(df['gene_center'] >= center - half) &
             (df['gene_center'] <= center + half)].copy()

    win['pip'] = win['pips_pop1'].fillna(0).clip(0, 1)
    # y-axis z: FUSION marginal .dat TWAS.Z if a dat dir is given (matched on
    # Ensembl for tissues / symbol for cell types), else FOCUS fine-mapping z.
    if fusion_dat_dir:
        chrom = int(win['chrom'].iloc[0])
        by_ens, by_sym = _fusion_z_map(fusion_dat_dir, chrom)

        def _fz(row):
            return by_ens.get((row['ensembl'], row['group']),
                              by_sym.get((row['symbol'], row['group']), np.nan))
        win['twas_z'] = win.apply(_fz, axis=1)
        miss = win['twas_z'].isna()
        if miss.any():
            print(f"[build] {int(miss.sum())} gene-models lacked a FUSION TWAS.Z "
                  f"-> fell back to FOCUS z")
            win.loc[miss, 'twas_z'] = win.loc[miss, 'twas_z_pop1']
        print(f"[build] z-source = FUSION .dat TWAS.Z ({fusion_dat_dir})")
    else:
        win['twas_z'] = win['twas_z_pop1']
        print("[build] z-source = FOCUS twas_z_pop1")
    win['in_cred_set'] = pd.to_numeric(win['in_cred_set_pop1'], errors='coerce').fillna(0).astype(int)
    # deterministic horizontal jitter so overlapping group-models are separable
    win['x_plot'] = (win['gene_center'] / 1e6
                     + np.random.RandomState(0).uniform(-0.006, 0.006, len(win)) * window_mb)

    win = win[CSV_COLS]
    win.to_csv(out_csv, index=False)
    print(f"[build] {block}\n[build] wrote {out_csv}: {win['ensembl'].nunique()} genes, "
          f"{win['group'].nunique()} groups, {len(win)} gene-models, "
          f"{int((win['in_cred_set']==1).sum())} credible models")
    return win


# ----------------------------------------------------------------------------
# stage 2: committed CSV -> figure
# ----------------------------------------------------------------------------
def _assign_rows(genes, window_mb, char_mb):
    """Greedy interval-graph row assignment so gene-track labels don't overlap."""
    rows, row_of = [], {}
    for _, r in genes.sort_values('tx_start').iterrows():
        w = max(len(r['symbol']) * char_mb, (r['g_stop'] - r['g_start']) / 1e6)
        c = (r['g_start'] + r['g_stop']) / 2 / 1e6
        lo, hi = c - w / 2, c + w / 2
        placed = False
        for i, last in enumerate(rows):
            if lo > last + 0.01 * window_mb:
                rows[i] = hi
                row_of[r['ensembl']] = i
                placed = True
                break
        if not placed:
            rows.append(hi)
            row_of[r['ensembl']] = len(rows) - 1
    return row_of, len(rows)


def regional_plot(csv_path, focal_ens, focal_name, panel_letter, group_label,
                  window_mb=1.0, out='panel.png', title='', extra_labels=None):
    win = pd.read_csv(csv_path)
    center = win.loc[win['ensembl'] == focal_ens, 'gene_center'].median()
    half = window_mb * 1e6 / 2

    genes = (win.groupby(['ensembl', 'symbol'])
                .agg(g_start=('tx_start', lambda s: s.iloc[0] + 1e6),
                     g_stop=('tx_stop', lambda s: s.iloc[0] - 1e6),
                     tx_start=('tx_start', 'first')).reset_index())

    # Figure height can be overridden (REGIONAL_FIG_H) to render a taller panel
    # for the combined Figure 3 bottom row, where a wide 12x7.4 panel would be
    # width-limited and look small. Standalone default stays 7.4.
    fig = plt.figure(figsize=(12, float(os.environ.get("REGIONAL_FIG_H", 7.4))))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.2, 1.6], hspace=0.06)
    ax = fig.add_subplot(gs[0])
    axg = fig.add_subplot(gs[1], sharex=ax)

    # ---- main scatter ----
    for z in (-1.96, 1.96):   # nominal two-sided p=0.05, matching the A/B Miami plots
        ax.axhline(z, ls='--', lw=0.8, color='grey', alpha=0.6)
    ax.axhline(0, lw=0.6, color='black', alpha=0.4)
    rest = win[win['ensembl'] != focal_ens]
    sc = ax.scatter(rest['x_plot'], rest['twas_z'], c=rest['pip'], cmap='viridis',
                    vmin=0, vmax=1, s=24 + 130 * rest['pip'], edgecolor='white',
                    linewidth=0.3, alpha=0.9, zorder=3)
    fo = win[win['ensembl'] == focal_ens]                       # red triangle = focal-gene identity
    ax.scatter(fo['x_plot'], fo['twas_z'], c=fo['pip'], cmap='viridis', vmin=0, vmax=1,
               marker='^', s=70 + 200 * fo['pip'], edgecolor='red', linewidth=1.7, zorder=5)
    cred = win[win['in_cred_set'] == 1]                          # black ring = credible set
    ax.scatter(cred['x_plot'], cred['twas_z'], s=70 + 180 * cred['pip'],
               facecolor='none', edgecolor='black', linewidth=1.6, zorder=6)
    cbar = fig.colorbar(sc, ax=[ax, axg], pad=0.01, fraction=0.045)
    cbar.set_label('FOCUS PIP', fontsize=12)

    # label strongly-supported credible models: "GENE (group)"
    lab = cred[cred['pip'] > 0.5].sort_values('x_plot')
    pos = lab[lab['twas_z'] >= 0].reset_index(drop=True)
    for i, r in pos.iterrows():
        left = (i == 0)
        ax.annotate(f"{r['symbol']}\n({r['group'].replace('_', ' ')})", (r['x_plot'], r['twas_z']),
                    xytext=(-22 if left else 26, 30 + 4 * i), textcoords='offset points',
                    ha='right' if left else 'left', va='bottom', fontsize=8.5, fontweight='bold',
                    arrowprops=dict(arrowstyle='-', lw=0.8, color='black'), zorder=8)
    for _, r in lab[lab['twas_z'] < 0].iterrows():
        ax.annotate(f"{r['symbol']}\n({r['group'].replace('_', ' ')})", (r['x_plot'], r['twas_z']),
                    xytext=(-55, 22), textcoords='offset points', ha='right', va='center',
                    fontsize=8.5, fontweight='bold',
                    arrowprops=dict(arrowstyle='-', lw=0.8, color='black'), zorder=8)

    # optional manual call-outs (grey italic) for specific gene x group dots,
    # e.g. a strong-Z / low-PIP point worth pointing out. format: (ensembl, group, text)
    for ens, grp, text in (extra_labels or []):
        m = win[(win['ensembl'] == ens) & (win['group'] == grp)]
        if not len(m):
            continue
        r = m.iloc[0]
        up = r['twas_z'] >= 0
        ax.annotate(text.replace('\\n', '\n'), (r['x_plot'], r['twas_z']),
                    xytext=(34, 22 if up else -34), textcoords='offset points',
                    ha='center', va='bottom' if up else 'top', fontsize=8.5,
                    fontstyle='italic', color='#333333',
                    arrowprops=dict(arrowstyle='-', lw=0.8, ls='--', color='#555555'), zorder=8)

    ax.set_ylabel('TWAS Z-score', fontsize=13)
    # Title is sized so that, after this (wide) panel is shrunk into the
    # combined Figure 3 grid (~0.34x), its title text matches the other panels.
    # The leading panel letter is drawn separately (below) as a big bold glyph so
    # it matches the prominent A/B/C/D letters instead of being title-sized.
    ax.set_title(title or f"{focal_name} locus — {group_label}",
                 fontsize=20, fontweight='bold', loc='left')
    # Big standalone panel letter, dropped onto the title's line (instead of
    # floating at the very top) so the panel has no wasted top band and sits
    # closer to the plot in the combined figure.
    fig.text(0.012, 0.915, panel_letter, fontsize=40, fontweight='bold',
             ha='left', va='top')
    ax.text(0.99, 0.03, f"{window_mb:g} Mb window (±{window_mb/2:g} Mb around {focal_name})\n"
            f"{win['ensembl'].nunique()} genes · {win['group'].nunique()} "
            f"{group_label.lower()} · {len(win)} gene-"
            f"{group_label.lower().rstrip('s').replace(' ', '-')} pairs into FOCUS",
            transform=ax.transAxes, ha='right', va='bottom', fontsize=9.5, color='dimgrey')
    ax.tick_params(labelbottom=False)
    ax.legend(handles=[
        Line2D([0], [0], marker='^', color='none', markerfacecolor='none',
               markeredgecolor='red', markeredgewidth=1.7, markersize=12,
               label=f'{focal_name} model (one per {group_label.lower().rstrip("s")})'),
        Line2D([0], [0], marker='o', color='none', markerfacecolor='none',
               markeredgecolor='black', markeredgewidth=1.6, markersize=12,
               label='In FOCUS 90% credible set')],
        loc='upper left', fontsize=9.5, frameon=True,
        title='outline = marker · fill color = PIP', title_fontsize=8.5)

    # ---- gene track ----
    char_mb = window_mb * 0.013
    row_of, nrows = _assign_rows(genes, window_mb, char_mb)
    for _, r in genes.iterrows():
        is_f = r['ensembl'] == focal_ens
        y = -row_of[r['ensembl']]
        s_mb, e_mb = r['g_start'] / 1e6, r['g_stop'] / 1e6
        axg.plot([s_mb, e_mb], [y, y], lw=4 if is_f else 2.4,
                 color='red' if is_f else '#444', solid_capstyle='round', zorder=3)
        axg.annotate(r['symbol'], ((s_mb + e_mb) / 2, y + 0.3), ha='center', va='bottom',
                     fontsize=9.5 if is_f else 8.5, color='red' if is_f else 'black',
                     fontweight='bold' if is_f else 'normal')
    axg.set_ylim(-nrows - 0.2, 1.2)
    axg.set_yticks([])
    for sp in ('left', 'right', 'top'):
        axg.spines[sp].set_visible(False)
    axg.set_xlabel(f"Chromosome {int(win['chrom'].iloc[0])} position (Mb)", fontsize=13)
    axg.set_xlim(center / 1e6 - half / 1e6, center / 1e6 + half / 1e6)

    fig.savefig(out, dpi=300, bbox_inches='tight')
    fig.savefig(os.path.splitext(out)[0] + '.pdf', bbox_inches='tight')
    print(f"[plot] saved {out} (+ .pdf)")


# ----------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--focus-dir', default='.', help='dir with results*/ and gene_annotation.txt.gz')
    p.add_argument('--focal-ens', required=True, help='focal gene Ensembl id (unversioned), e.g. ENSG00000073605')
    p.add_argument('--focal-name', required=True, help='focal gene symbol, e.g. GSDMB')
    p.add_argument('--panel-letter', default='C')
    p.add_argument('--group-label', default='Tissues', help='"Tissues" or "Cell types"')
    p.add_argument('--results-glob', default='results/*/*.focus.tsv',
                   help="glob (relative to --focus-dir) of raw FOCUS block .tsv; "
                        "use 'results_ct/*/*.focus.tsv' for cell types")
    p.add_argument('--annotation', default=None,
                   help='gene_annotation.txt.gz (default: <focus-dir>/data/gene_annotation.txt.gz '
                        'if present, else <focus-dir>/gene_annotation.txt.gz)')
    p.add_argument('--window-mb', type=float, default=1.0)
    p.add_argument('--fusion-dat-dir', default=None,
                   help='dir of FUSION .dat files (twas_<ctx>_chr<N>.dat or '
                        'TWAS_metaGWAS_<ctx>_chr<N>.dat). If set, dot y = FUSION '
                        'marginal TWAS.Z instead of FOCUS twas_z (matched on '
                        'Ensembl for tissues / symbol for cell types).')
    p.add_argument('--csv', default='panel_C_data.csv', help='committed dataset path')
    p.add_argument('--out', default='panel_C_17q21_tissues.png')
    p.add_argument('--plot-only', action='store_true', help='skip rebuild; plot from existing --csv')
    p.add_argument('--title', default='')
    p.add_argument('--annotate-point', action='append', default=[], metavar='ENS:GROUP:TEXT',
                   help="grey italic call-out on a specific gene x group dot; repeatable. "
                        r"use \n in TEXT for line breaks. "
                        r"e.g. 'ENSG00000073605:esophagus_mucosa:GSDMB (esophagus mucosa)\nstrongest Z, not credible'")
    a = p.parse_args()
    extra = [tuple(s.split(':', 2)) for s in a.annotate_point]

    csv = a.csv if os.path.isabs(a.csv) else (
        a.csv if a.focus_dir == '.' else os.path.join(a.focus_dir, a.csv))
    if not a.plot_only:
        ann = a.annotation
        if ann is None:
            cand = os.path.join(a.focus_dir, 'data', 'gene_annotation.txt.gz')
            ann = cand if os.path.exists(cand) else os.path.join(a.focus_dir, 'gene_annotation.txt.gz')
        build_dataset(a.focus_dir, a.focal_ens, a.results_glob, a.window_mb, csv, ann,
                      fusion_dat_dir=a.fusion_dat_dir)
    regional_plot(csv, a.focal_ens, a.focal_name, a.panel_letter, a.group_label,
                  window_mb=a.window_mb, out=a.out, title=a.title, extra_labels=extra)


if __name__ == '__main__':
    main()
