#!/bin/bash

# In this job, we're running GRIEVOUS realign to harmonize variant IDs and
# allele orientations across two different cohorts. Within our study, one
# example of this would be the GACRS and 1000 Genomes Phase 3 cohorts

################################################################################
#SBATCH --job-name="grievous_realign"
#SBATCH --chdir="/path/to/logs/grievous_realign"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=1-22
#SBATCH --export=ALL
# --------------------------------------------------------------------------- #

echo "==================== SLURM JOB DEBUGGING INFORMATION ===================="
echo "[ SLURM Submit Directory (SLURM_SUBMIT_DIR) ]"
echo "$SLURM_SUBMIT_DIR"
echo "---------------------- SLURM Basic Job Information ----------------------"
echo "SLURM Job Name            [SLURM_JOB_NAME]      : ${SLURM_JOB_NAME}"
echo "SLURM Job Identifier      [SLURM_JOB_ID]        : ${SLURM_JOB_ID}"
echo "++ SLURM Array Identifier [SLURM_ARRAY_JOB_ID]  : ${SLURM_ARRAY_JOB_ID}"
echo "++ SLURM Array Task ID    [SLURM_ARRAY_TASK_ID] : ${SLURM_ARRAY_TASK_ID}"
echo "SLURM Job Partition       [SLURM_JOB_PARTITION] : ${SLURM_JOB_PARTITION}"
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
#                analysis; GRIEVOUS DB is named with 1KG first 
#                (flip order for QC database convention)
# ============================================================================ #
# COHORT INFORMATION for GTEx (Expression Reference)
#   Study      : Genotype-Tissue Expression (GTEx)
#   Note       : Used as reference eQTL panel; aligned to GACRS as primary
# ============================================================================ #
# COHORT INFORMATION for OneK1K (Single-Cell eQTL Reference)
#   Study      : OneK1K
#   Note       : Used as single-cell eQTL reference; included as all samples
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
GENOME_ASSEMBLY=hg38
MAPPING_FILE=${SLURM_SUBMIT_DIR}/grievous_pvar-mapping.txt

# GRIEVOUS database name for the pairwise alignment
GRIEVOUS_DB_NAME_COHORT1_COHORT2="cohort1-merged-cohort2"

echo ""
echo "Genome Assembly: ${GENOME_ASSEMBLY}"
echo "====================="
echo "Generating the Following GRIEVOUS Databases:"
echo "+ COHORT 1 merged COHORT 2 : ${GRIEVOUS_DB_NAME_COHORT1_COHORT2}"
echo ""

echo "-------------------------------------------------------------------------"
echo "Initializing CHROMOSOME Parameters for Each Dataset"
echo "-------------------------------------------------------------------------"
CHROMOSOME=$SLURM_ARRAY_TASK_ID

INPUT_FN_COHORT1=$( find $INPUT_DIR_COHORT1/*.pvar | grep "chr${CHROMOSOME}\." )
INPUT_FN_COHORT2=$( find $INPUT_DIR_COHORT2/*.pvar | grep "chr${CHROMOSOME}\." )

echo "Evaluating Chromosome: chr${CHROMOSOME}"
echo "============================"
echo "+ Input COHORT1 File : ${INPUT_FN_COHORT1}"
echo "+ Input COHORT2 File : ${INPUT_FN_COHORT2}"
echo ""

echo "-------------------------------------------------------------------------"
echo "Running GRIEVOUS REALIGN for COHORT1 + COHORT2"
echo "-------------------------------------------------------------------------"

mkdir -p ${OUTPUTS_DIR}/cohort1-merged-cohort2/cohort1
mkdir -p ${OUTPUTS_DIR}/cohort1-merged-cohort2/cohort2

pixi run grievous realign \
    --file "${INPUT_FN_COHORT1}" \
    --assembly "${GENOME_ASSEMBLY}" \
    --database "${GRIEVOUS_DB_NAME_COHORT1_COHORT2}" \
    -w "${OUTPUTS_DIR}/cohort1-merged-cohort2/cohort1" \
    --mapping "${MAPPING_FILE}" \
    --verbose --return_all -c "##"

pixi run grievous realign \
    --file "${INPUT_FN_COHORT2}" \
    --assembly "${GENOME_ASSEMBLY}" \
    --database "${GRIEVOUS_DB_NAME_COHORT1_COHORT2}" \
    -w "${OUTPUTS_DIR}/cohort1-merged-cohort2/cohort2" \
    --mapping "${MAPPING_FILE}" \
    --verbose --return_all -c "##"

echo "Done Running GRIEVOUS REALIGN for COHORT 1 + COHORT 2 !"

echo "-------------------------------------------------------------------------"
echo "Done with script !"
