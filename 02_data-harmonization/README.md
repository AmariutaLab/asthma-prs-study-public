# Data harmonization and pre-processing

This directory contains the code and associated notes for conducting data harmonization and pre-processing steps for the pediatric asthma polygenic risk score (PRS) project.

The pipeline harmonizes genotype data across two TOPMed whole-genome sequencing (WGS) cohorts along with several reference datasets: the Genetics of Asthma in Costa Rica Study (GACRS; phs001365, Freeze 10b); the Childhood Asthma Management Program (CAMP; phs001726, Freeze 10b); the 1000 Genomes Project Phase 3 Cohort; the Genotype-Tissue Expression (GTEx) Cohort; and the OneK1K Cohort. This general process produces a single merged, variant-aligned dataset for downstream analysis between the two cohorts of choice to merge.

## Tools used

Every tool invoked by the scripts in this directory is listed below so you can
decide which to install. Versions in parentheses are what the study used.

| Tool | Used for | Where |
|------|----------|-------|
| [PLINK 1.9](https://www.cog-genomics.org/plink/) (`plink`) | ID/sex/parent/phenotype updates, sample filtering, cross-cohort `--bmerge`, `--pca` | Steps 01, 07, 08 |
| [PLINK2](https://www.cog-genomics.org/plink/2.0/) (`plink2`) | VCF export, PGEN conversion + variant-ID standardization, `--rm-dup`, `--alt1-allele`, `--update-name`, `--extract`, `--pmerge-list` | Steps 02, 03, 05, 07 |
| [bcftools](https://samtools.github.io/bcftools/) (v1.21) | `view` (sample subset / biallelic filter), `norm` (split multiallelics), `+fixref` plugin (allele orientation vs FASTA), `sort`, `concat`, `index` | Steps 02, 03 |
| [CrossMap](https://crossmap.sourceforge.net/) (`CrossMap vcf`) | hg19→hg38 liftover — **OneK1K only** | Step 02b |
| [GRIEVOUS](https://github.com/biona001/Ghostknife) (`grievous realign`/`merge`/`intersect`) | pairwise variant-ID and allele harmonization across cohorts | Steps 04, 06 |
| `tar`, `rsync` | extract dbGaP `*.tar.gz`, stage VCFs | Step 01 |
| `awk`, `grep`, `cut`, `sort`, `find` (coreutils) | build sample/SNP lists, split intersection files by chromosome, merge-list assembly | Steps 01–07 |
| [SLURM](https://slurm.schedmd.com/) (*optional*) | per-chromosome array jobs (`--array=1-22`) | Steps 01–08 |

Environment/package manager used in the study: [pixi](https://pixi.sh/) (drives
the GRIEVOUS + CrossMap environment); [conda](https://docs.conda.io/) works too.
`bcftools` is on bioconda; the `+fixref` plugin ships with bcftools ≥1.10.

## Environment setup

Two environments — genotype tools, and GRIEVOUS/CrossMap (which need `numpy<2`).

```bash
# genotype tools
conda create -n asthma-demo -c bioconda -c conda-forge \
    plink=1.90b7.7 plink2 bcftools=1.21
# PLINK 1.9 v1.9.0-b.7.7, PLINK 2.0 v2.0.0-a.6.9, bcftools 1.21 (ships +fixref)

# GRIEVOUS + CrossMap (grievous imports fail under NumPy 2.x — pin numpy<2)
conda create -n grievous -c conda-forge python=3.10 'numpy<2' pandas
conda activate grievous
pip install grievous CrossMap        # grievous 0.1.5, CrossMap 0.7.0
```

The GRCh38 FASTA (`bcftools +fixref`) and `hg19ToHg38.over.chain.gz` (CrossMap)
are large public downloads, not packages — see **Data requirements** below.

## Data requirements

| Data | Status | Source / notes |
|------|--------|----------------|
| GACRS WGS (`phs001365`, Freeze 10b) | **Controlled-access** | dbGaP application required |
| CAMP WGS (`phs001726`, Freeze 10b) | **Controlled-access** | dbGaP application required |
| GTEx / OneK1K genotypes | **Controlled-access** | dbGaP / EGA application required |
| 1000 Genomes Project Phase 3 | **Public** | reference panel; used as the demo genotype source |
| `GRCh38_full_analysis_set_plus_decoy_hla.fa` (+`.fai`) | **Public, large (~3 GB)** | indexed hg38 FASTA (e.g. 1000G / UCSC); required by `bcftools +fixref` — not committed |
| `hg19ToHg38.over.chain.gz` | **Public** | UCSC chain file; required by `CrossMap` (OneK1K only) — not committed |

No individual-level study data is committed to this repository.

## Reproducible demo (smoke test)

> [!NOTE]
> The demo is **planned** for this module and not yet committed (see the repo's
> build status). When added, it will live in `02_realign-variants/sample_data/`
> with a `run_demo.sh` runner, following the same pattern as
> [`03_ancestry-analysis/`](../03_ancestry-analysis/).

Design: Step 01 (extract dbGaP archives) and Steps 02–03 (`bcftools +fixref`
against the multi-GB GRCh38 FASTA) depend on controlled-access data and an
un-committable reference, so the demo will **start from pre-built PLINK2 PGEN
files** and exercise the harmonization core — **Steps 04–08** (`grievous
realign` → `plink2` reorientation → `grievous merge`/`intersect` → cross-cohort
`--bmerge` → `--pca`). Two pseudo-cohorts are derived from public 1000 Genomes
chr22 (split into disjoint sample sets, with a subset of alleles deliberately
flipped so GRIEVOUS has real work to do). All demo inputs will be labeled as a
**real 1000G subset**; nothing individual-level from the study cohorts is used.

> [!IMPORTANT]
> GRIEVOUS installed from PyPI currently imports NumPy 1.x-compiled extensions
> and fails under NumPy 2.x. Pin `numpy<2` in the demo environment
> (`conda install 'numpy<2'` or `pip install 'numpy<2'`).

## General Analysis Steps

### Step 01 — Extract dbGaP Archive and Filter Samples

The first script, `01_plink_untar-grab-samples.sh`, processes raw `*.tar.gz` archives downloaded from dbGaP for GACRS and CAMP, filtering each dataset down to the samples of interest and adding the necessary sample metadata. This script runs as a single batch job that loops over chromosomes 1–22.

#### Script Overview

1. **Untar the dbGaP archive** — Extract the raw `*.tar.gz` archive to scratch storage using `tar`, then `rsync` the decompressed VCF files to our Lustre storage on EXPANSE for easier access
2. **Update sample IDs** — Run `plink --update-ids` to assign correct IID/FID values to each sample using a pre-prepared mapping file
3. **Add sample metadata** — Update sex (`--update-sex`), parental information (`--update-parents`), and phenotype (`--make-pheno`) for each sample; require that all retained samples have sex information (`--must-have-sex`)
4. **Filter to samples of interest** — Retain only cases and unrelated controls (`--keep`) and exclude duplicate variants (`--exclude *.dupvar`)

The per-chromosome reference files (ID maps, sex/parent/phenotype files, and the cases-and-controls sample list) are expected to be located in `REFERENCE_DIR` and named `{COHORT}_update-plink-{ids,sex,parents,phenotype}.txt` and `{COHORT}_cases-and-controls.txt`.

The final output is a set of per-chromosome PLINK binary filesets (`*.updated.cases-unrel-controls-only`) written to `OUTPUTS_DIR`.

---

### Step 02 — `bcftools` QC and Reference Allele Fixing

Two variants of this step exist depending on the genome build of the input data. Both run as SLURM array jobs (one task per chromosome, `--array=1-22`).

#### Step 02a — `bcftools` QC for GRCh38-formatted Files

The script `02a_bcftools_clean-chromosomes.sbatch` applies bcftools quality control to the PLINK binary filesets produced in Step 01, normalizing multiallelic sites and fixing allele orientations against the GRCh38 reference.

##### Script Overview

1. **Recode to VCF** — Convert the PLINK binary fileset to VCF format using `plink2 --export vcf`
2. **Extract samples** — Subset the VCF to the samples of interest using `bcftools view --samples-file`
3. **Decompose multiallelic sites** — Split multiallelic sites into biallelic records using `bcftools norm --multiallelics -any`, tagging the original records with `OLD_MULTIALLELIC_VAR`
4. **Fix allele orientations** — Align REF/ALT alleles against the GRCh38 reference FASTA using `bcftools +fixref --mode flip --discard`
    * Flips REF/ALT columns and genotypes for non-ambiguous SNPs
    * Discards sites that cannot be resolved against the reference
5. **Sort and index** — Sort the output using `bcftools sort` and write an index with `--write-index`

The final output is a set of per-chromosome sorted, indexed BCF files (`*.sorted.bcf`) written to `OUTPUTS_DIR`.

#### Step 02b — CrossMap Liftover + `bcftools` QC for GRCh37-formatted Files

The script `02b_bcftools_clean-chromosomes-with-crossmap.sbatch` follows the same process as Step 02a, but inserts a CrossMap liftover step between sample extraction and bcftools normalization to convert coordinates from hg19 to hg38.

##### Script Overview

1. **Recode to VCF** and **extract samples** (same as Step 02a)
2. **Liftover from hg19 to hg38** — Run `CrossMap vcf` with the `hg19ToHg38.over.chain.gz` chain file and the GRCh38 reference FASTA
3. **Decompose multiallelic sites**, **fix allele orientations**, and **sort and index** (same as Step 02a)

The final output is a set of per-chromosome sorted, indexed BCF files (`*.sorted.bcf`) written to `OUTPUTS_DIR`.

---

### Step 03 — Convert to PLINK2 PGEN Format

The script `03_plink_convert-chromosomes.sbatch` takes the per-chromosome BCF files from Step 02 and produces per-chromosome PLINK2 PGEN files for input into the GRIEVOUS variant harmonization pipeline.

#### Script Overview

1. **Concatenate per-chromosome BCF files** into a single merged BCF using `bcftools concat --allow-overlaps`, then sort the result
2. **Filter for biallelic SNPs only** using `bcftools view --min-alleles 2 --max-alleles 2 --types snps`
3. **Split by chromosome and convert to PLINK2 PGEN format** using `plink2 --make-pgen`, standardizing variant IDs to `CHR:POS:REF:ALT` format with `--set-all-var-ids @:#:$r:$a`

The final output is a set of per-chromosome PLINK2 PGEN files with normalized, biallelic SNPs and standardized variant IDs written to `OUTPUTS_DIR`.

---

### Step 04 — GRIEVOUS Realign

The script `02_realign-variants/01_grievous_realign-variants.sbatch` runs `grievous realign` on each pairwise combination of cohorts to harmonize variant IDs and allele orientations in preparation for merging. This script runs as a SLURM array job (one task per chromosome, `--array=1-22`).

For each chromosome, `grievous realign` is run once per cohort within each pairwise combination, producing per-chromosome reorientation files that specify how alleles should be flipped to achieve cross-cohort consistency.

> [!NOTE]
> A `grievous_pvar-mapping.txt` file is required to define custom column mappings for the `.pvar` input files. This file is included in the `02_realign-variants/` directory.

The final output is a set of per-chromosome `Reorientation/` directories, each containing allele flip instructions (`NoDuplicates_ReorientRefAlleleThisWay_CHR{N}.tsv`) and variant ID remapping files (`NoDuplicates_ReorientIndex_CHR{N}.tsv`).

---

### Step 05 — Apply GRIEVOUS Allele Reorientation

The script `02_realign-variants/02_plink2_harmonization-grievous-realigned-variants.sbatch` applies the reorientation output from Step 04 to the actual PLINK2 genotype files. This script runs as a SLURM array job (one task per chromosome, `--array=1-22`).

#### Script Overview

1. **Remove duplicate variants** — Run `plink2 --rm-dup force-first` on both cohort files to ensure each variant position appears only once
2. **Reorient alleles** — Apply the GRIEVOUS allele flip instructions using `plink2 --alt1-allele 'force'`
    * The `force` modifier prevents PLINK2 from erroring when a "known" reference allele is changed
3. **Update variant IDs** — Remap variant identifiers to the harmonized GRIEVOUS index using `plink2 --update-name`

The final output is a set of per-chromosome realigned and reindexed PLINK2 PGEN files (`*.realigned.reindexed`) written to `OUTPUTS_DIR`.

---

### Step 06 — GRIEVOUS Merge and Intersect

The script `02_realign-variants/03_grievous_merge-intersect-variants.sbatch` aggregates the per-chromosome GRIEVOUS realign results and identifies the set of SNPs shared across each cohort pair.

#### Script Overview

1. **Merge per-chromosome realign outputs** — Run `grievous merge -a` on each cohort's realigned output directory to combine the per-chromosome results into a single merged variant list
2. **Identify overlapping SNPs** — Run `grievous intersect` on each cohort pair to produce a single SNP list containing only variants present in both cohorts

The final output for each cohort pair is an `AllIntersectingSNPs_{cohort1}-merged-{cohort2}.tsv` file containing the set of shared variant IDs written to `OUTPUTS_DIR`.

---

### Step 07 — Extract Intersecting SNPs and Merge Cohorts

The script `02_realign-variants/04_plink2_intersect-grievous-snps.sbatch` uses the intersection SNP lists from Step 06 to extract overlapping variants from each cohort's PLINK2 files per chromosome, then merges across chromosomes and across cohorts into a final merged dataset.

#### Script Overview

1. **Split the intersection SNP list by chromosome** — Partition the GRIEVOUS intersection file into per-chromosome SNP lists using `grep`
2. **Extract overlapping SNPs per chromosome** — Run `plink2 --extract` on each cohort's per-chromosome PLINK2 files using the chromosome-specific SNP list
3. **Merge across chromosomes within each cohort** — Run `plink2 --pmerge-list` to concatenate the per-chromosome PLINK2 files for each cohort into a single genome-wide file in both BED and PGEN format
4. **Merge across cohorts** — Use `plink --bmerge --real-ref-alleles` to merge the two cohorts' BED files into a single combined dataset, then convert back to PLINK2 PGEN format with `plink2 --make-pgen`

> [!NOTE]
> PLINK 1.9 is required for the cross-cohort merge step because PLINK2 does not currently support non-concatenation merges (i.e., merging datasets with different sample sets).

The final output is a merged PLINK BED + PGEN fileset containing all samples from both cohorts, restricted to the GRIEVOUS-intersecting SNPs.

---

### Step 08 — PCA on Combined Cohort

Finally, the script `02_realign-variants/05_plink_combined-cohorts-pca.sbatch` runs principal component analysis (PCA) on the final merged dataset to characterize population structure across all cohorts.

#### Script Overview

1. **Run PCA** — Execute `plink --pca 10 'header' 'tabs' 'var-wts'` on the merged PLINK binary fileset, computing 10 principal components with both sample eigenvectors and per-variant loadings

The final output is a set of PCA files (`*.eigenvec`, `*.eigenval`, `*.eigenvec.var`) written to `OUTPUTS_DIR` for use in downstream population stratification analysis.
