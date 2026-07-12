"""
LD reference heatmap (r^2) at the GSDMB / 17q21 asthma locus.

Renders a lower-triangle r^2 heatmap over the top-N cis-SNPs (ranked by max
|FUSION eQTL weight| across GTEx tissues) using the 1000 Genomes EUR
reference panel bundled with FUSION (TWAS/LDREF/1000G.EUR.17.{bed,bim,fam}).

Pipeline:
    1. Load per-tissue GSDMB weights CSV produced by
       ../../05_twas-analysis/figure_scripts/extract_gtex_weights.R
       -> columns: tissue, snp, weight
    2. Rank SNPs by max |weight| across tissues; take top N.
    3. Call plink twice:
         plink --bfile <ldref> --extract snp_list --r square
         plink --bfile <ldref> --extract snp_list --write-snplist
       (plink emits SNPs in bim / genomic order, which is also the axis order
        in the heatmap.)
    4. r^2 = plink r squared; mask upper triangle; render with matplotlib.

Outputs (see CONFIG):
    FOCUS/ld_reference_1kg.png     (300 dpi, transparent)
    FOCUS/ld_reference_1kg.pdf     (vector, editable text)
"""
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Editable text in SVG / PDF -- keep glyph data, don't outline to paths
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype']  = 42


# =============================================================================
# CONFIG -- edit these paths for your environment (local laptop vs server)
# =============================================================================
ROOT = Path(os.environ.get('HARTWELL_ROOT',
                           '/Users/nancyh/Desktop/hartwell'))

# Long CSV from extract_gtex_weights.R (tissue, snp, weight)
WEIGHTS_CSV = Path(os.environ.get(
    'GSDMB_WEIGHTS_CSV',
    Path(tempfile.gettempdir()) / 'gsdmb_blup4.csv'))

# FUSION 1000G LD reference prefix (chr17 for the GSDMB locus)
LDREF_BFILE = Path(os.environ.get(
    'LDREF_BFILE', ROOT / 'TWAS' / 'LDREF' / '1000G.EUR.17'))

# plink 1.9 binary
PLINK_BIN = Path(os.environ.get('PLINK_BIN', ROOT / 'plink'))

# Output directory
OUTDIR = Path(os.environ.get('FIG_OUTDIR', ROOT / 'FOCUS'))
OUTDIR.mkdir(parents=True, exist_ok=True)

# How many top-|weight| SNPs to include in the heatmap
TOP_N = int(os.environ.get('LD_TOP_N', '26'))

# Locus label for the plot title
LOCUS_GENE   = os.environ.get('LOCUS_GENE',   'GSDMB')
LOCUS_REGION = os.environ.get('LOCUS_REGION', '17q21')


# =============================================================================
# Step 1 -- pick top-N SNPs by max |weight| across tissues
# =============================================================================
if not WEIGHTS_CSV.exists():
    raise FileNotFoundError(
        f"Weights CSV not found: {WEIGHTS_CSV}\n"
        f"Produce it first with:\n"
        f"  Rscript 05_twas-analysis/figure_scripts/extract_gtex_weights.R"
    )

d = pd.read_csv(WEIGHTS_CSV)
if not {'tissue', 'snp', 'weight'} <= set(d.columns):
    raise ValueError(f"Expected columns tissue/snp/weight in {WEIGHTS_CSV}; "
                     f"got {list(d.columns)}")

# Rank SNPs by max |weight| across tissues, take top N
score = (d.assign(abs_w=d['weight'].abs())
          .groupby('snp')['abs_w'].max()
          .sort_values(ascending=False))
top_snps = score.head(TOP_N).index.tolist()
print(f"Selected top {len(top_snps)} SNPs from {WEIGHTS_CSV.name} "
      f"(pool: {d['snp'].nunique():,} unique SNPs)")


# =============================================================================
# Step 2 -- plink extract + LD matrix
# =============================================================================
tmpdir = Path(tempfile.mkdtemp(prefix='ld_ref_'))
snp_list = tmpdir / 'top_snps.txt'
snp_list.write_text('\n'.join(top_snps) + '\n')

ld_prefix   = tmpdir / 'ld'
match_prefix = tmpdir / 'match'


def _run_plink(*extra_args):
    cmd = [str(PLINK_BIN), '--bfile', str(LDREF_BFILE),
           '--extract', str(snp_list), '--silent', *extra_args]
    print('  $ ' + ' '.join(cmd))
    subprocess.run(cmd, check=True)


print(f"Running plink over 1000G reference {LDREF_BFILE}.{{bed,bim,fam}} ...")
_run_plink('--r', 'square', '--out', str(ld_prefix))
_run_plink('--write-snplist', '--out', str(match_prefix))

ld_path    = ld_prefix.with_suffix('.ld')
match_path = match_prefix.with_suffix('.snplist')
if not ld_path.exists() or not match_path.exists():
    raise RuntimeError(f"plink did not produce {ld_path} or {match_path} -- "
                       f"check its log at {ld_prefix}.log / {match_prefix}.log")


# =============================================================================
# Step 3 -- load + plot lower-triangle r^2
# =============================================================================
r = np.loadtxt(ld_path)
matched_snps = [line.strip() for line in open(match_path) if line.strip()]

n_dropped = len(top_snps) - len(matched_snps)
if n_dropped > 0:
    print(f"  {n_dropped} SNP(s) not in the 1000G panel; "
          f"proceeding with {len(matched_snps)} matched SNPs.")

r2 = r ** 2
mask = np.triu(np.ones_like(r2, dtype=bool), k=1)   # hide upper triangle
A = np.ma.masked_where(mask, r2)

fig, ax = plt.subplots(figsize=(5.2, 5.2), dpi=300)
cmap = plt.cm.Reds.copy()
cmap.set_bad('white')
im = ax.imshow(A, cmap=cmap, vmin=0, vmax=1)
ax.set_xticks(range(len(matched_snps)))
ax.set_yticks(range(len(matched_snps)))
ax.set_xticklabels(matched_snps, rotation=90, fontsize=5)
ax.set_yticklabels(matched_snps, fontsize=5)
ax.set_title(
    f'LD reference (1000G EUR)\n'
    fr'$\it{{{LOCUS_GENE}}}$ locus $\cdot$ {LOCUS_REGION} $\cdot$ $r^2$',
    fontsize=10, fontweight='bold')
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label(r'$r^2$', fontsize=9)
cb.ax.tick_params(labelsize=7)
plt.tight_layout()

out_png = OUTDIR / 'ld_reference_1kg.png'
out_pdf = OUTDIR / 'ld_reference_1kg.pdf'
plt.savefig(out_png, bbox_inches='tight', dpi=300)
plt.savefig(out_pdf, bbox_inches='tight')

mean_offdiag = r2[~np.eye(len(matched_snps), dtype=bool)].mean()
print()
print(f'Saved: {out_png}')
print(f'Saved: {out_pdf}')
print(f'mean off-diagonal r^2: {mean_offdiag:.3f}')
