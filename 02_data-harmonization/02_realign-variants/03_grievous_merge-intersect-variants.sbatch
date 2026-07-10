#!/bin/bash

# In this job, we're running GRIEVOUS MERGE and GRIEVOUS INTERSECT to identify
# overlapping SNPs across all cohort pairings (i.e., our cohorts gacrs or camp
# versus reference datasets such as 1KG, gtex, and onek1k) in preparation for
# the downstream PLINK2-based variant filtering

################################################################################
#SBATCH --job-name="grievous_merge-intersect-variants"
#SBATCH --chdir="/path/to/logs/grievous_merge-intersect-variants"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
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

INPUT_DIR=${PROJECT_DIR}/path/to/user/output/from/script/02
OUTPUT_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p $OUTPUT_DIR

################################################################################
echo "-------------------------------------------------------------------------"
echo "Initializing Project Environment with PIXI"
echo "-------------------------------------------------------------------------"
cd $PROJECT_DIR
eval "$(pixi shell-hook --environment data-preprocessing)"
pixi info
echo ""

echo "-------------------------------------------------------------------------"
echo "Running GRIEVOUS MERGE for COHORT 1-COHORT 2 Datasets"
echo "-------------------------------------------------------------------------"

MERGE_COHORT1_COHORT2=${INPUT_DIR}/cohort1-merged-cohort2
pixi run grievous merge -a ${MERGE_COHORT1_COHORT2}/cohort1
pixi run grievous merge -a ${MERGE_COHORT1_COHORT2}/cohort2
echo "+ Completed GRIEVOUS MERGE for COHORT 1 and COHORT 2 !"

echo "-------------------------------------------------------------------------"
echo "Running GRIEVOUS INTERSECT for COHORT 1-COHORT 2 Datasets"
echo "-------------------------------------------------------------------------"

pixi run grievous intersect \
    --write_path ${OUTPUT_DIR} \
    --cohorts ${MERGE_COHORT1_COHORT2}/cohort1 ${MERGE_COHORT1_COHORT2}/cohort2 \
    -n AllIntersectingSNPs_cohort1-merged-cohort2

echo "-------------------------------------------------------------------------"
echo "Done with script !"
echo ""
