#!/bin/bash

# In this job, we're converting our files to PLINK2 format for input
# into the GRIEVOUS pipeline (i.e., requires PLINK2 over PLINK1.9)

################################################################################
#SBATCH --job-name="plink2_convert-chromosomes"
#SBATCH --chdir="/path/to/logs/plink2_convert-chromosomes"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --export=ALL
# --------------------------------------------------------------------------- #

# Output debugging information to the output logs for diagnosing bugs later ...
JOB_SCRATCH=/scratch/$USER/job_$SLURM_JOB_ID

export BCFTOOLS_PLUGINS=${BCFTOOLS_PLUGINS_DIR}
echo "==================== SLURM JOB DEBUGGING INFORMATION ===================="
echo "[ SLURM Submit Directory (SLURM_SUBMIT_DIR) ]"
echo "$SLURM_SUBMIT_DIR"
echo "---------------------- SLURM Basic Job Information ----------------------"
echo "SLURM Job Name            [SLURM_JOB_NAME]      : ${SLURM_JOB_NAME}"
echo "SLURM Job Identifier      [SLURM_JOB_ID]        : ${SLURM_JOB_ID}"
echo "++ SLURM Array Identifier [SLURM_ARRAY_JOB_ID]  : ${SLURM_ARRAY_JOB_ID}"
echo "++ SLURM Array Task ID    [SLURM_ARRAY_TASK_ID] : ${SLURM_ARRAY_TASK_ID}"
echo "SLURM Job Partition       [SLURM_JOB_PARTITION] : ${SLURM_JOB_PARTITION}"
echo "SLURM Job Scratch Location                      : ${JOB_SCRATCH}"
echo "-------------------- SLURM Job Resources Information --------------------"
echo "SLURM Job Number of Nodes          [SLURM_JOB_NUM_NODES]   : ${SLURM_JOB_NUM_NODES}"
echo "SLURM Job Total MEM per Node       [SLURM_MEM_PER_NODE]    : ${SLURM_MEM_PER_NODE}"
echo "SLURM Job Number of Tasks per Node [SLURM_NTASKS_PER_NODE] : ${SLURM_NTASKS_PER_NODE}"
echo "SLURM Job Number of CPUs           [SLURM_NPROCS]          : ${SLURM_NPROCS}"
echo "========================================================================="
echo ""

# End script automatically with errors (-e) and return error code (-o pipefail)
set -e -o pipefail
################################################################################

# ============================================================================ #
# COHORT INFORMATION for GACRS Freeze 10b
#   Study      : Genetics of Asthma in Costa Rica Study (GACRS)
#   dbGaP Acc. : phs001365
#   Dataset    : phg002069.v1 — TOPMed WGS Asthma CostaRica v6 frz10 SNP
#   Access Tier: DS-ASTHMA-IRB-MDS-RD
#   Samples    : Children ages 6–14 with asthma (cases + unrelated controls)
# ============================================================================ #
# COHORT INFORMATION for CAMP Freeze 10b
#   Study      : Childhood Asthma Management Program (CAMP)
#   dbGaP Acc. : phs001726
#   Dataset    : phg002060.v1 — TOPMed WGS CAMP v3 frz10 SNP
#   Access Tier: DS-AST-COPD
#   Samples    : Children ages 5–12 with asthma (cases + unrelated controls)
# ============================================================================ #

PROJECT_DIR=/path/to/user/project
BCFTOOLS_PLUGINS_DIR=/path/to/user/bcftools-1.21/plugins

COHORT_NAME=gacrs # or camp, gtex, onek1k
COHORT_NAME_UC=$( echo $COHORT_NAME | tr '[:lower:]' '[:upper:]' )

# Define the glob pattern that we can use to match the per-chromosome BCF files
# from the previous step, relative to the sorted/ subdirectory where our input
# files are located to use with the find UNIX utility
BCF_GLOB_PATTERN="*freeze.10b*.bcf"

# Define a shorter cohort label to use as the prefix for final output files
OUTPUT_PREFIX=gacrs_frz10b
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output/from/script/02
INPUT_DIR=${OUTPUTS_DIR}/sorted # inputs are in a subdirectory from the previous script

################################################################################

TMP_MERGED_UNSORTED="${JOB_SCRATCH}"/tmp_${OUTPUT_PREFIX}_allchr.updated.cases-unrel-controls-only.unsorted.bcf
TMP_MERGED_SORTED="${JOB_SCRATCH}"/tmp_${OUTPUT_PREFIX}_allchr.updated.cases-unrel-controls-only.sorted.bcf
MERGED_NORMED="${OUTPUTS_DIR}"/${OUTPUT_PREFIX}_allchr.updated.cases-unrel-controls-only.normed.bcf

echo "Gathering Chromosome Files for ${COHORT_NAME_UC} Dataset to Merge ..."
find "${INPUT_DIR}"/"${BCF_GLOB_PATTERN}" | sort -V > "${JOB_SCRATCH}"/tmp_mergelist.txt

echo "Merging Chromosome Files into Single Entity for ${COHORT_NAME_UC} Dataset ..."
${PROJECT_DIR}/./bcftools concat --allow-overlaps \
    --file-list "${JOB_SCRATCH}"/tmp_mergelist.txt \
    --output "${TMP_MERGED_UNSORTED}"
${PROJECT_DIR}/./bcftools sort \
    "${TMP_MERGED_UNSORTED}" \
    --output "${TMP_MERGED_SORTED}" \
    --write-index
echo "Done Merging Chromosome Files into Single Entity for ${COHORT_NAME_UC} Dataset !"

echo "Filtering for Biallelic SNPs ONLY for ${COHORT_NAME_UC} Dataset ..."
${PROJECT_DIR}/./bcftools view \
    --min-alleles 2 \
    --max-alleles 2 \
    --types snps \
    "${TMP_MERGED_SORTED}" \
    --output "${MERGED_NORMED}"
${PROJECT_DIR}/./bcftools index "${MERGED_NORMED}"
echo "Done Filtering for Biallelic SNPs ONLY for ${COHORT_NAME_UC} Dataset !"

echo "Filtering on Chromosome by Chromosome Basis ..."
for CHROMOSOME in {1..22}; do
    echo "+ Extracting chr${CHROMOSOME} as PLINK2-Formatted File ..."
    ${PROJECT_DIR}/./plink2 \
        --bcf "${MERGED_NORMED}" \
        --chr "${CHROMOSOME}" \
        --snps-only just-acgt --max-alleles 2 \
        --set-all-var-ids @:#:\$r:\$a \
        --new-id-max-allele-len 1000 \
        --make-pgen \
        --out "${OUTPUTS_DIR}"/${OUTPUT_PREFIX}_chr"${CHROMOSOME}".updated.cases-unrel-controls-only.normed
done

# Cleanup
rm "${JOB_SCRATCH}"/tmp_mergelist.txt
rm "${TMP_MERGED_UNSORTED}"*
rm "${TMP_MERGED_SORTED}"*

echo "Done Filtering on Chromosome by Chromosome Basis !"
echo ""
echo "-------------------------------------------------------------------------"
echo ""
echo "Done with script !"
echo ""
