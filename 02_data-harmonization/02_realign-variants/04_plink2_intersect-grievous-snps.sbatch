#!/bin/bash

# In this job, we're using the SNP intersection lists from GRIEVOUS INTERSECT
# to extract overlapping variants from each cohort's PLINK2 files per
# chromosome, then merging across chromosomes and across cohort pairs into our
# final PLINK1.9 and PLINK2 datasets for downstream analysis

################################################################################
#SBATCH --job-name="plink2_intersect-grievous-snps"
#SBATCH --chdir="/path/to/logs/plink2_intersect-grievous-snps"
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
# COHORT INFORMATION for GACRS Freeze 10b (Primary Cohort)
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
# COHORT INFORMATION for 1000 Genomes Project (1KG; Ancestry Reference)
#   Study      : 1000 Genomes Project (Phase 3)
#   Note       : Used as an ancestry reference panel for population structure
#                analysis; harmonized to GACRS as the primary cohort
# ============================================================================ #
# COHORT INFORMATION for GTEx (Expression Reference)
#   Study      : Genotype-Tissue Expression (GTEx)
#   Note       : Used as reference eQTL panel; harmonized to GACRS as primary
# ============================================================================ #
# COHORT INFORMATION for OneK1K (Single-Cell eQTL Reference)
#   Study      : OneK1K
#   Note       : Used as single-cell eQTL reference; processed as all samples
#                and as a subset filtered for autoimmune-relevant cell types
# ============================================================================ #

################################################################################
# Note: For these variables in the script, a lot of the paths and names here
#       will be different based on your environment setup! We've removed a lot
#       of the lab-specific path information, but the code logic is the same ...
################################################################################

PROJECT_DIR=/path/to/user/project

# ---------------------------------------------------------------------------- #
# Handle the input files first!
INTERSECT_DIR=${PROJECT_DIR}/path/to/user/output/from/script/03

INTERSECT_FN_COHORT1_COHORT2=${INTERSECT_DIR}/AllIntersectingSNPs_cohort1-merged-cohort2.tsv
INPUT_DIR_COHORT1=${PROJECT_DIR}/path/to/user/output/from/script/02 # previously realigned
INPUT_DIR_COHORT2=${PROJECT_DIR}/path/to/plink2/pgen/cohort2

echo "SNP Intersection Files from GRIEVOUS INTERSECT:"
echo "-----------------------------------------------"
echo "+ Intersection File for COHORT 1 + COHORT 2 : ${INTERSECT_FN_COHORT1_COHORT2}"
echo ""

echo "Input PLINK File Locations for Downstream Filtering:"
echo "----------------------------------------------------"
echo "+ COHORT 1 : ${INPUT_DIR_COHORT1}"
echo "+ COHORT 2 : ${INPUT_DIR_COHORT2}"
echo ""

# ---------------------------------------------------------------------------- #
# Handle the output files next!
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p ${OUTPUTS_DIR}/cohort1-merged-cohort2

# Define the intermediate output filenames, then the final merged output
FINAL_OUTPUT_PLINK_COHORT1_COHORT2=${OUTPUTS_DIR}/cohort1-merged-cohort2/cohort2.grievous-filt.cohort1-merged-cohort2
FINAL_OUTPUT_PLINK_COHORT1=${OUTPUTS_DIR}/cohort1-merged-cohort2/cohort1.grievous-filt.cohort1-merged-cohort2
FINAL_OUTPUT_PLINK_MERGED_COHORT1_COHORT2=${OUTPUTS_DIR}/cohort1.merged.cohort2.grievous-filt

################################################################################
TMP_OUTPUT_DIR_COHORT1_COHORT2="${JOB_SCRATCH}"/tmp-outdir_cohort1-merged-cohort2.cohort2-filt
TMP_OUTPUT_DIR_COHORT1="${JOB_SCRATCH}"/tmp-outdir_cohort1-merged-cohort2.cohort1-filt

mkdir -p "$TMP_OUTPUT_DIR_COHORT1_COHORT2"
mkdir -p "$TMP_OUTPUT_DIR_COHORT1"

echo "========================================================================="
echo ""
sort -V ${INTERSECT_FN_COHORT1_COHORT2} > "${JOB_SCRATCH}"/tmp-intersect_cohort1-merged-cohort2.txt
echo "Sorted SNP Intersection Files for COHORT 1 + COHORT 2 !"
echo ""

echo "========================================================================="
echo ""
echo "Processing All Chromosomes for COHORT 1 + COHORT 2 ..."
echo ""

for CHROMOSOME in {1..22}; do
    echo "+ Extracting SNPs for chr${CHROMOSOME} ..."
    INPUT_PLINK_COHORT1=$( find "${INPUT_DIR_COHORT1}"/*.pvar | grep "chr${CHROMOSOME}\." | sed 's/.pvar//' )
    INPUT_PLINK_COHORT2=$( find "${INPUT_DIR_COHORT2}"/*.pvar | grep "chr${CHROMOSOME}\." | sed 's/.pvar//' )

    echo "---------------------------------------------------------------------"
    echo "+ >>> Input COHORT 1 PLINK File : ${INPUT_PLINK_COHORT1}"
    echo "+ >>> Input COHORT 2 PLINK File : ${INPUT_PLINK_COHORT2}"
    echo "---------------------------------------------------------------------"
    echo ""

    echo "+ Splitting SNP Intersection File for chr${CHROMOSOME} ..."
    OVERLAPPING_COHORT1_COHORT2="${JOB_SCRATCH}"/tmp-overlapping-snps.cohort1-merged-cohort2.chr"${CHROMOSOME}".txt
    grep "^${CHROMOSOME}:" "${JOB_SCRATCH}"/tmp-intersect_cohort1-merged-cohort2.txt > "${OVERLAPPING_COHORT1_COHORT2}"

    echo "---------------------------------------------------------------------"
    echo "+ >>> Separated Overlapping SNPs in COHORT 1 + COHORT 2 for chr${CHROMOSOME} : N=$( wc -l "${OVERLAPPING_COHORT1_COHORT2}" ) variants"
    echo "---------------------------------------------------------------------"
    echo ""

    echo "+ Extracting SNPs with PLINK2 for COHORT 1 + COHORT 2 Files ..."
    OUTPUT_PLINK_COHORT1_COHORT2=$( basename "$INPUT_PLINK_COHORT2" )
    OUTPUT_PLINK_COHORT1_COHORT2=${OUTPUT_PLINK_COHORT1_COHORT2//normed/grievous-filt}
    ${PROJECT_DIR}/./plink2 --pfile "$INPUT_PLINK_COHORT2" \
        --extract "$OVERLAPPING_COHORT1_COHORT2" \
        --make-pgen --out "${TMP_OUTPUT_DIR_COHORT1_COHORT2}/${OUTPUT_PLINK_COHORT1_COHORT2}"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    OUTPUT_PLINK_COHORT1=$( basename "$INPUT_PLINK_COHORT1" )
    OUTPUT_PLINK_COHORT1=${OUTPUT_PLINK_COHORT1//normed.realigned.reindexed/grievous-filt}
    ${PROJECT_DIR}/./plink2 --pfile "$INPUT_PLINK_COHORT1" \
        --extract "$OVERLAPPING_COHORT1_COHORT2" \
        --make-pgen --out "${TMP_OUTPUT_DIR_COHORT1}/${OUTPUT_PLINK_COHORT1}"

    echo "---------------------------------------------------------------------"
    echo ""
done

echo "========================================================================="
echo "Generating MERGE LIST for COHORT 1-COHORT 2 Files"
echo ""

find "${TMP_OUTPUT_DIR_COHORT1_COHORT2}"/*.pvar | sed 's/.pvar//' > "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort2.txt
find "${TMP_OUTPUT_DIR_COHORT1}"/*.pvar | sed 's/.pvar//' > "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort1.txt

echo "Merge List for COHORT 1-COHORT 2 Files"
echo "========================================"
echo "COHORT 1 Files"
echo "----------------------------------------"
cat "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort1.txt
echo ""
echo "----------------------------------------"
echo "COHORT 2 Files"
echo "----------------------------------------"
cat "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort2.txt
echo ""

echo "========================================================================="
echo "Merging PLINK2 Files Across Chromosomes for COHORT 1-COHORT 2 Files"
echo ""

# We are generating plink1.9 files as well because plink2 can't handle
# non-concatenation merges at the minute (i.e., merges between different sample sets)

${PROJECT_DIR}/./plink2 \
    --pmerge-list "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort2.txt \
    --make-bed --out "$FINAL_OUTPUT_PLINK_COHORT1_COHORT2"
${PROJECT_DIR}/./plink2 \
    --pmerge-list "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort2.txt \
    --make-pgen --out "$FINAL_OUTPUT_PLINK_COHORT1_COHORT2"

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo ""

${PROJECT_DIR}/./plink2 \
    --pmerge-list "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort1.txt \
    --make-bed --out "$FINAL_OUTPUT_PLINK_COHORT1"
${PROJECT_DIR}/./plink2 \
    --pmerge-list "${JOB_SCRATCH}"/tmp-merge.cohort1-merged-cohort2.cohort1.txt \
    --make-pgen --out "$FINAL_OUTPUT_PLINK_COHORT1"

echo "------------------------------------------------------------"
echo "+ Completed PLINK2 Merge for COHORT 1 + COHORT 2s !"
echo "------------------------------------------------------------"
echo ""

echo "========================================================================="
echo "Merging PLINK1.9 Files Across Cohorts for COHORT 1-COHORT 2 Files"

${PROJECT_DIR}/./plink --bfile "$FINAL_OUTPUT_PLINK_COHORT1_COHORT2" \
    --bmerge "$FINAL_OUTPUT_PLINK_COHORT1" --real-ref-alleles \
    --make-bed --out "$FINAL_OUTPUT_PLINK_MERGED_COHORT1_COHORT2"

echo "========================================================================="
echo "Converting Merged PLINK1.9 Files to PLINK2 Format for COHORT 1-COHORT 2 Files"

${PROJECT_DIR}/./plink2 --bfile "$FINAL_OUTPUT_PLINK_MERGED_COHORT1_COHORT2" \
    --make-pgen --out "$FINAL_OUTPUT_PLINK_MERGED_COHORT1_COHORT2"

echo "========================================================================="
echo "Done with script !"
echo ""
