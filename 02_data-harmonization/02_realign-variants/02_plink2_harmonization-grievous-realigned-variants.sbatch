#!/bin/bash

# In this job, we're applying GRIEVOUS realign output to harmonize allele
# orientations across the reference datasets (e.g., 1KG, GTEx, OneK1K) aligned
# to our cohorts (e.g., GACRS or CAMP)

################################################################################
#SBATCH --job-name="plink2_grievous-harmonization"
#SBATCH --chdir="/path/to/logs/plink2_grievous-harmonization"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1-22
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

# Input directories for each dataset (i.e., in PLINK2 PGEN format)
INPUT_DIR_COHORT1=${PROJECT_DIR}/path/to/plink2/pgen/cohort1
INPUT_DIR_COHORT2=${PROJECT_DIR}/path/to/plink2/pgen/cohort2

# Specify the GRIEVOUS realign output directories from the previous script
# as they contain the per-chromosome Reorientation/ subdirectories
GRIEVOUS_REALIGN_DIR=${PROJECT_DIR}/path/to/user/output/from/script/01
REALIGN_DIR_COHORT1_COHORT2=${GRIEVOUS_REALIGN_DIR}/cohort1-merged-cohort2/cohort1/Reorientation

OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p $OUTPUTS_DIR

################################################################################
echo "-------------------------------------------------------------------------"
echo "Initializing Project Environment with PIXI"
echo "-------------------------------------------------------------------------"
cd $PROJECT_DIR
eval "$(pixi shell-hook --environment data-preprocessing)"
pixi info
echo ""

echo "-------------------------------------------------------------------------"
echo "Initializing GRIEVOUS Parameters"
echo "-------------------------------------------------------------------------"
echo ""
echo "GRIEVOUS REALIGN Output Directories"
echo "==================================="
echo "+ COHORT 1-COHORT 2 Database : ${REALIGN_DIR_COHORT1_COHORT2}"
echo ""

echo "-------------------------------------------------------------------------"
echo "Initializing CHROMOSOME Parameters for Each Dataset"
echo "-------------------------------------------------------------------------"
CHROMOSOME=$SLURM_ARRAY_TASK_ID

INPUT_FN_COHORT1=$( find $INPUT_DIR_COHORT1/*.pvar | grep "chr${CHROMOSOME}\." )
INPUT_FN_COHORT2=$( find $INPUT_DIR_COHORT2/*.pvar | grep "chr${CHROMOSOME}\." )

# GRIEVOUS output filenames follow a fixed convention (CHR-suffixed TSV files)
REALIGN_ALLELE_FN=NoDuplicates_ReorientRefAlleleThisWay_CHR${CHROMOSOME}.tsv
REALIGN_INDEX_FN=NoDuplicates_ReorientIndex_CHR${CHROMOSOME}.tsv

echo "Evaluating Chromosome: chr${CHROMOSOME}"
echo "============================"
echo "+ Input COHORT1 File : ${INPUT_FN_COHORT1}"
echo "+ Input COHORT2 File : ${INPUT_FN_COHORT2}"
echo ""
echo "- GRIEVOUS Allele Flip File       : ${REALIGN_ALLELE_FN}"
echo "- GRIEVOUS Allele Flip Index File : ${REALIGN_INDEX_FN}"
echo ""

echo "-------------------------------------------------------------------------"
echo "Removing DUPLICATE SNPs from Unaligned PLINK Files for chr${CHROMOSOME}"
echo "-------------------------------------------------------------------------"

BASENAME_COHORT1=$( basename "$INPUT_FN_COHORT1" )
BASENAME_COHORT2=$( basename "$INPUT_FN_COHORT2" )

mkdir -p "${JOB_SCRATCH}"/cohort1_tmp
mkdir -p "${JOB_SCRATCH}"/cohort2_tmp

${PROJECT_DIR}/./plink2 \
    --pfile "${INPUT_FN_COHORT1}" \
    --rm-dup force-first \
    --make-pgen \
    --out "${JOB_SCRATCH}"/cohort1_tmp/"${BASENAME_COHORT1}".nodup
echo "+ Completed Removal of DUPLICATE SNPs for COHORT 1 Dataset !"

${PROJECT_DIR}/./plink2 \
    --pfile "${INPUT_FN_COHORT2}" \
    --rm-dup force-first \
    --make-pgen \
    --out "${JOB_SCRATCH}"/cohort1_tmp/"${BASENAME_COHORT2}".nodup
echo "+ Completed Removal of DUPLICATE SNPs for COHORT 2 Dataset !"

echo "-------------------------------------------------------------------------"
echo "Reorienting FLIPPED SNPs Identified with GRIEVOUS REALIGN for chr${CHROMOSOME}"
echo "-------------------------------------------------------------------------"

# --alt1-allele 'force' prevents plink2 from erroring on "known" reference allele
# changes; see: https://www.cog-genomics.org/plink/2.0/data#ref_allele

mkdir -p ${OUTPUTS_DIR}/cohort1-merged-cohort2

${PROJECT_DIR}/./plink2 \
    --pfile "${JOB_SCRATCH}"/cohort1_tmp/"${BASENAME_COHORT1}".nodup \
    --alt1-allele 'force' "${REALIGN_DIR_COHORT1_COHORT2}/${REALIGN_ALLELE_FN}" \
    --make-pgen \
    --out "${JOB_SCRATCH}"/cohort1_tmp/"${BASENAME_COHORT1}".realigned
${PROJECT_DIR}/./plink2 \
    --pfile "${JOB_SCRATCH}"/cohort1_tmp/"${BASENAME_COHORT1}".realigned \
    --max-alleles 2 \
    --update-name "${REALIGN_DIR_COHORT1_COHORT2}/${REALIGN_INDEX_FN}" \
    --make-pgen \
    --out "${OUTPUTS_DIR}"/cohort1-merged-cohort2/"${BASENAME_COHORT1}".realigned.reindexed
echo "+ Completed Reorientation for COHORT 1 Dataset against COHORT 2 Dataset !"

echo "-------------------------------------------------------------------------"
echo ""
echo "Done with script !"
echo ""
