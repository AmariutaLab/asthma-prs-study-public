# Polygenic risk score (PRS) construction with PRS-CS and PRS-CSx

This folder contains the code and associated notes for constructing the
**polygenic risk score (PRS)** with the [`PRS-CS`](https://github.com/getian107/PRScs)
and [`PRS-CSx`](https://github.com/getian107/PRScsx) toolkits for the pediatric
asthma polygenic risk score (PRS) project.

`PRS-CS` places continuous-shrinkage priors on SNP effect sizes using a single
(European-ancestry) LD reference panel; `PRS-CSx` is the cross-population
extension that couples multiple ancestry panels in a single model. Both consume
the meta-analysis GWAS summary statistics produced in
[`01_meta-analysis/`](../01_meta-analysis/) and produce posterior SNP effect
sizes that are then applied to a validation genotype target with `plink --score`.

## Tools used

Every tool invoked by the scripts in this directory is listed below so you can
decide which to install. Versions in parentheses are what the study used.

| Tool | Used for | Where |
|------|----------|-------|
| [PRS-CS](https://github.com/getian107/PRScs) (`PRScs.py`) | posterior SNP effect sizes, single EUR LD panel, per chromosome, across φ (`1e-2`, `1e-4`, `auto`) | `01_prscs-scripts/01` |
| [PRS-CSx](https://github.com/getian107/PRScsx) (`PRScsx.py`) | cross-population posterior effect sizes over multiple LD panels | `02_prscsx-scripts/01` |
| Python 3 + [`scipy`](https://scipy.org/), [`numpy`](https://numpy.org/), [`h5py`](https://www.h5py.org/) | required by PRS-CS / PRS-CSx (h5py reads the `ldblk_*.hdf5` panels) | `*_scripts/01` |
| [PLINK 1.9](https://www.cog-genomics.org/plink/) (`plink`) | `--merge-list` across chromosomes, `--score <file> 2 4 6 sum center` to compute PRS | `*_scripts/03` |
| `cat`, `find`, `sort`, `sed`, `awk` (coreutils) | merge per-chromosome posterior files, assemble merge lists, reformat summary statistics | `*_scripts/02`, `03` |
| [SLURM](https://slurm.schedmd.com/) (*optional*) | array jobs over (chromosome × φ) parameter grid — see `prscs_run-parameters.txt` / `prscsx_run-parameters.txt` | `*_scripts/01`, `03` |

Environment/package manager used in the study: [pixi](https://pixi.sh/) (the
`prscs-construction` environment); [conda](https://docs.conda.io/) works too.
PRS-CS / PRS-CSx are cloned from GitHub (not pip packages):

```bash
git clone https://github.com/getian107/PRScs.git
git clone https://github.com/getian107/PRScsx.git
conda create -n asthma-demo -c conda-forge python=3.11 scipy numpy h5py
conda install -n asthma-demo -c bioconda plink
```

## Environment setup

PRS-CS / PRS-CSx are GitHub clones (not pip packages); the Python env needs a
matched `numpy`/`h5py` (a clean env avoids ABI errors).

```bash
# python env — h5py must match numpy, so pin numpy<2 in a fresh env
conda create -n prscs-demo -c conda-forge python=3.10 'numpy<2' scipy h5py pandas
conda activate prscs-demo

# scoring
conda install -c bioconda plink=1.90b7.7          # PLINK 1.9 v1.9.0-b.7.7

# toolkits (version 1.1.0)
git clone https://github.com/getian107/PRScs
git clone https://github.com/getian107/PRScsx
```

## Data requirements

| Data | Status | Source / notes |
|------|--------|----------------|
| Meta-analysis GWAS summary statistics (`SNP A1 A2 BETA SE`) + total sample size | **Derived** | output of [`01_meta-analysis/`](../01_meta-analysis/); an example is in `01_meta-analysis/sample_data/expected_outputs/` |
| Validation genotype target (`.bed/.bim/.fam`, per chromosome) | **Controlled-access** for the real study cohorts | GACRS/CAMP dbGaP; 1000 Genomes used as the public demo target |
| LD reference panel `ldblk_1kg_eur/` (+`snpinfo_1kg_hm3`) | **Public, large (~4.5 GB)** | [PRS-CS panels](https://github.com/getian107/PRScs#getting-started); not committed |
| Additional-ancestry LD panels (`ldblk_1kg_eas`, `_afr`, …) | **Public, large** | required only for PRS-CSx cross-population runs |

No individual-level study data is committed to this repository.

## General Analysis Steps

Both `01_prscs-scripts/` and `02_prscsx-scripts/` follow the same three-stage
shape; PRS-CSx differs by taking **one summary-statistics file and one LD panel
per ancestry** and adding a `--pop` argument.

### Step 01 — Run PRS-CS / PRS-CSx

`01_prscs_run-meta-analysis-with-cohort.sbatch` (and its PRS-CSx twin) runs the
tool once per `(chromosome, φ)` combination, looked up from the accompanying
`*_run-parameters.txt` array file. Inputs:

- `--sst_file` — meta-analysis summary statistics formatted as `SNP A1 A2 BETA SE`
- `--n_gwas` — total GWAS sample size (cases + controls)
- `--bim_prefix` — the validation target's PLINK `.bim` prefix (defines the SNP set to score)
- `--ref_dir` — the `ldblk_*` LD reference panel
- `--phi` — global shrinkage parameter (fixed `1e-2` / `1e-4`, or `auto`)

Output per `(chromosome, φ)`: a posterior effect-size table
`..._pst_eff_a1_b0.5_phi{φ}_chr{N}.txt` (columns: `CHR SNP BP A1 A2 BETA`).

### Step 02 — Merge posterior effect sizes across chromosomes

`02_bash_merge-prscs-results.sh` concatenates the per-chromosome posterior files
into one genome-wide file per φ (`..._phi{φ}_allchr.txt`).

### Step 03 — Score the validation target

`03_plink_generate-scores-with-prscs-result.sbatch` merges the per-chromosome
target PLINK files (`plink --merge-list`) and applies each φ's posterior weights
with `plink --score <file> 2 4 6 sum center` (column 2 = SNP, 4 = A1, 6 = BETA),
yielding one `.profile` PRS per φ for downstream evaluation.

## Reproducible demo (smoke test)

A committed fixture + runner exercise **PRS-CS and PRS-CSx end-to-end** on public
data:

```bash
cd 04_prscs-prscsx-construction
PRSCS_DIR=/path/to/PRScs PRSCSX_DIR=/path/to/PRScsx \
  PYTHON=python PLINK=plink ./run_demo.sh
```

It runs Steps 01–03 on **chromosome 22, 250 HapMap3 SNPs**: PRS-CS + PRS-CSx MCMC
→ merge → `plink --score`, using a **real 1000G-EUR LD reference** rebuilt from
public chr22 genotypes, a **real 1000G-EUR target** (100 samples), and a **FAKE,
simulated** GWAS sumstats file — all in [`sample_data/`](sample_data/) (~310 KB).
MCMC is seed-fixed, so the posteriors reproduce the committed reference exactly
(max \|Δbeta\| = 0). See [`sample_data/README.md`](sample_data/README.md).

What is real vs. simulated:

- **Summary statistics** — **FAKE / simulated** `SNP A1 A2 BETA SE`, clearly
  labeled; the *format* matches `01_meta-analysis` (in production that real
  meta-GWAS is the input).
- **Validation target** — **REAL** public 1000 Genomes Phase 3 EUR chr22 subset.
- **LD reference** — a tiny chr22-only `ldblk` HDF5 **rebuilt from the public
  1000G demo genotypes** in the exact PRS-CS format, so the tools run end-to-end.
  The official ~4.5 GB `ldblk_1kg_eur` panel is documented (link above) for real runs.

> [!NOTE]
> Requires a clone of [`getian107/PRScs`](https://github.com/getian107/PRScs) and
> [`getian107/PRScsx`](https://github.com/getian107/PRScsx), PLINK 1.9
> (v1.9.0-b.7.7), and Python 3 with `scipy`/`numpy`/`h5py` (use a clean env so
> `h5py` matches `numpy`, e.g. `conda create -c conda-forge python=3.10 'numpy<2'
> scipy h5py`).

The canonical `.sbatch` / `.sh` scripts remain the reference HPC versions;
`run_demo.sh` drops the SLURM / lab paths and points at the demo inputs.
