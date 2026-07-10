#!/bin/bash

# This script filters the inputted plink file to HapMap3 variants (rsID format)
# obtained from the Alkes Group. You can obtain this file in compressed format 
# from: https://zenodo.org/records/7773502

# We also perform some quality control filtering on the plink files; we are
# inputting files that have already been pre-combined with 1000 Genomes samples
# (e.g., GACRS+1KG, etc.)

################################################################################
# Note: For these variables in the script, a lot of the paths and names here
#       will be different based on your environment setup! We've removed a lot
#       of the lab-specific path information, but the code logic is the same ...
################################################################################

PROJECT_DIR=/path/to/user/project
JOB_SCRATCH=/scratch/$USER/job_$SLURM_JOB_ID # directory for temporary files on EXPANSE

HAPMAP3_VARIANTS=${PROJECT_DIR}/data/w_hm3.snplist

COHORT_NAME=gacrs # or camp, gtex, onek1k
INPUT_PLINK=${PROJECT_DIR}/data/cohort-data/name-of-plink-file-to-process # .bed/.bim/.fam

OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p "$OUTPUTS_DIR"

################################################################################
# Skip the header in the w_hm3.snplist and save to a temporary file
TMP_HAPMAP3_RSIDS=${JOB_SCRATCH}/tmp_rsids.txt
awk 'NR > 1 {print $1}' "${HAPMAP3_VARIANTS}" > "${TMP_HAPMAP3_RSIDS}"

################################################################################
# Note : We had some metadata in another directory that we wanted to add to our
#        files (e.g., phenotype information, individual sex). This is optional!
################################################################################

BASENAME=$( basename "${INPUT_PLINK}" )

# Filter for HapMap3 variants from the plink file
"${PROJECT_DIR}"/./plink \
    --bfile "${INPUT_PLINK}" \
    --update-sex "${JOB_SCRATCH}"/${COHORT_NAME}_update-plink-sex.txt \
    --make-pheno "${JOB_SCRATCH}"/${COHORT_NAME}_update-plink-phenotype.txt '*' \
    --extract "${TMP_HAPMAP3_RSIDS}" \
    --allow-no-sex \
    --make-bed --out "${JOB_SCRATCH}"/"${BASENAME}"

# Perform QC for variant missingness (--geno 0.05) and MAF (--maf 0.01)
"${PROJECT_DIR}"/./plink \
    --bfile "${JOB_SCRATCH}"/"${BASENAME}" \
    --geno 0.05 --maf 0.01 \
    --make-bed --out "${JOB_SCRATCH}"/"${BASENAME}".variant-qc

# Perform LD pruning (based on --indep-pairwise 50 10 0.1)
"${PROJECT_DIR}"/./plink \
    --bfile "${JOB_SCRATCH}"/"${BASENAME}".variant-qc \
    --indep-pairwise 50 10 0.1 \
    --out "${JOB_SCRATCH}"/"${BASENAME}".ld-qc
"${PROJECT_DIR}"/./plink \
    --bfile "${JOB_SCRATCH}"/"${BASENAME}".variant-qc \
    --extract "${JOB_SCRATCH}"/"${BASENAME}".ld-qc.prune.in \
    --make-bed --out "${OUTPUTS_DIR}"/"${BASENAME}".filt-hm3

# Perform cleanup of temporary files from our scratch directory
rm "${JOB_SCRATCH}"/"${BASENAME}".variant-qc*
rm "${JOB_SCRATCH}"/"${BASENAME}".ld-qc*
rm "${JOB_SCRATCH}"/"${BASENAME}"* # any files we missed

################################################################################
echo ""
echo "Done with script !"
echo ""
