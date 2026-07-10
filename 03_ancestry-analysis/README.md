# Ancestry analysis with ADMIXTURE

This directory contains the code and associated notes for performing ancestry analysis with `ADMIXTURE` for the pediatric asthma polygenic risk score (PRS) and polygenic transcriptomic risk score (PTRS) project.

## Tools used

Every tool invoked by the scripts in this directory is listed below so you can
decide which to install. Versions in parentheses are what the study used.

| Tool | Used for | Where |
|------|----------|-------|
| [PLINK 1.9](https://www.cog-genomics.org/plink/) (`plink`, v1.9.0-b.7.7) | HapMap3 `--extract`, variant QC (`--geno`, `--maf`), LD pruning (`--indep-pairwise`), optional `--update-sex`/`--make-pheno` | Step 01 |
| [ADMIXTURE](https://dalexander.github.io/admixture/) (v1.3.0) | unsupervised ancestry estimation with cross-validation (`--cv`, `-j`) | Step 02 |
| `awk` (coreutils) | strip the header off `w_hm3.snplist` | Step 01 |
| [SLURM](https://slurm.schedmd.com/) (*optional*) | array-parallel execution across *K* values | Step 02 |

Package managers used in the study (either works): [pixi](https://pixi.sh/) or
[conda](https://docs.conda.io/). PLINK and ADMIXTURE are both on bioconda:
`conda create -n asthma-demo -c bioconda plink admixture`.

## Environment setup

Both tools are on bioconda:

```bash
conda create -n asthma-demo -c bioconda -c conda-forge \
    plink=1.90b7.7 admixture=1.3.0
conda activate asthma-demo
# PLINK 1.9 v1.9.0-b.7.7, ADMIXTURE 1.3.0
```

## Data requirements

| Data | Status | Source / notes |
|------|--------|----------------|
| Genotypes merged with the 1000 Genomes reference (`.bed/.bim/.fam`) | **Controlled-access** for the real study cohorts (GACRS/CAMP dbGaP) | produced by [`02_data-harmonization/`](../02_data-harmonization/) |
| `w_hm3.snplist` (HapMap3 SNP list, rsID format) | **Public** | Alkes Group, [Zenodo 7773502](https://zenodo.org/records/7773502) |
| Per-sample sex / phenotype metadata (*optional*) | Cohort-specific | only needed for the optional `--update-sex` / `--make-pheno` flags |

No individual-level study data is committed to this repository.

## Reproducible demo (smoke test)

A self-contained demonstration that runs both steps end-to-end on public data
lives alongside the canonical scripts:

```bash
cd 03_ancestry-analysis
PLINK=plink ADMIXTURE=admixture ./run_demo.sh     # binaries from PATH by default
```

`run_demo.sh` mirrors the logic of `01_plink_filter-hm3-variants.sh` and
`02_admixture_calculate-ancestry-unsupervised.sbatch` but drops the SLURM / lab
paths and points at the committed demo inputs in [`sample_data/`](sample_data/):

- `sample_data/demo_1kg.chr22.{bed,bim,fam}` — **REAL** 1000 Genomes Project
  Phase 3 EUR genotypes, chr22 only, downsampled to 100 samples / 6,000 SNPs.
- `sample_data/w_hm3_chr22_subset.snplist` — **REAL** HapMap3 list (as above),
  subset to the chr22 rsIDs in the demo.

> [!NOTE]
> The demo genotypes are single-ancestry (1000G EUR), so this is a smoke test of
> the *pipeline*, not a biological result — cross-validation error will favour
> *K*=1. A genuinely multi-ancestry input produces meaningful structure from the
> same commands. See [`sample_data/README.md`](sample_data/README.md) for exact
> provenance and regeneration steps.

## General Analysis Steps

### Step 01 — Filter Genotype Data to HapMap3 Variants

The first script, `01_plink_filter-hm3-variants.sh`, filters a PLINK binary fileset (i.e., `.bed/.bim/.fam`) down to HapMap3 variants and applies quality control filters prior to ancestry analysis (Step #2). Input files are expected to already be merged with 1000 Genomes samples (i.e., Genetics of Asthma in Costa Rica (GACRS) Study individuals + 1000 Genomes Project individuals).

#### Script Overview

1. **Extract HapMap3 Variant IDs** — Generate a temporary copy of the `w_hm3.snplist` file without the header for input to `plink`

> [!NOTE]
> If you do not already have the HapMap3 SNP list, the Alkes Group has [uploaded it to Zenodo](https://zenodo.org/records/7773502).

2. **Filter to HapMap3 Variants** — Run PLINK on the input files with `--extract` to retain only HapMap3 rsIDs
    * (*Optional*) During this step, we also updated per-sample sex (`--update-sex`) and phenotype (`--make-pheno`) information from metadata, but you can ignore these flags if this information is already integrated with your data or if it is not necessary
3. **Perform Variant Quality Control** — Apply filters for genotype missingness (`--geno 0.05`) and minor allele frequency (`--maf 0.01`)
    * We exclude variants with a genotype missingness rate of >5% to reduce bias from non-random missingness (which can correlate with batch, ancestry, or genotype itself), ensure reliable ancestry inference, and avoid artificial discordance across merged cohorts where variants may be well-genotyped in one dataset but not another
    * We exclude variants with a MAF <1% as rare variants contribute little information to ancestry inference (i.e., ADMIXTURE relies on allele frequency differences between populations, and rare variants observed too infrequently across samples provide unstable frequency estimates that can distort cluster assignment)
4. **Conduct LD Pruning** — Prune variants in linkage disequilibrium (LD) using a sliding window (`--indep-pairwise 50 10 0.1`; window = 50 SNPs, step = 10 SNPs, r2 threshold = 0.1) where the pruned-in variant list is used to produce the final filtered fileset
    * We remove correlated variants to ensure that remaining SNPs are approximately independent (i.e., ADMIXTURE assumes independence between variants and retaining correlated SNPs may bias cluster assignments)

The final output is a PLINK fileset named `<input-basename>.filt-hm3` written to `OUTPUTS_DIR`.

---

### Step 02 — Unsupervised ADMIXTURE Analysis

Finally, the script `02_admixture_calculate-ancestry-unsupervised.sbatch` runs ADMIXTURE in unsupervised mode across a range of *K* values (i.e., number of ancestral populations) on the HapMap3-filtered, LD-pruned PLINK fileset produced in Step #1. We ran this script with the SLURM workload manager to enable parallel cross-validation across multiple *K* values in a single job (i.e., each array task here evaluates a single *K*)

#### Script Overview

1. **Set the Value of *K* to Evaluate** — The `SLURM_ARRAY_TASK_ID` is used directly as the *K* value, so array indices map 1:1 to *K* values (e.g., array `1-5` runs *K*=1 through *K*=5)
2. **Run ADMIXTURE with Cross-Validation** — Execute ADMIXTURE with `--cv` (cross-validation) and `-j` (parallelism set to `SLURM_NTASKS_PER_NODE`)

> [!NOTE]
> ADMIXTURE writes output files to the current working directory, so the script also `cd`'s into `OUTPUTS_DIR` before running

Per *K*, we generate the following output files:

- `${BASENAME}.${K}.Q` — The Ancestry Proportion Matrix (i.e., one row per sample)
- `${BASENAME}.${K}.P` — Allele Frequency Matrix (i.e., one row per variant)
- Cross-Validation Error printed to `stdout` (i.e., used to select optimal *K*)
