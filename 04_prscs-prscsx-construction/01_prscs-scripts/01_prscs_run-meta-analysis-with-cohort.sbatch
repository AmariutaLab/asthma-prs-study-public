#!/bin/bash

################################################################################
#SBATCH --job-name="prscs_run-meta-analysis-with-cohort"
#SBATCH --chdir="/path/to/logs/prscs_run-meta-analysis-with-cohort"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --array=1-66
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=36G
#SBATCH --time=08:00:00
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
echo "Starting Time: $( date )"
echo ""

# End script automatically with errors (-e) and return error code (-o pipefail)
set -e -o pipefail

################################################################################
# Note: For these variables in the script, a lot of the paths and names here
#       will be different based on your environment setup! We've removed a lot
#       of the lab-specific path information, but the code logic is the same ...
################################################################################

PROJECT_DIR=/path/to/user/project
ARRAY_CONFIG_FILE=${SLURM_SUBMIT_DIR}/prscs_run-parameters.txt

# We're using the PLINK file to be used as a validation target dataset
# containing the SNPs of interest within our summary statistics reference
INPUTS_DIR=${PROJECT_DIR}/data/cohort-data/name-of-plink-file-to-process
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p ${OUTPUTS_DIR}

################################################################################
echo "Initializing Environment with PIXI"
cd ${PROJECT_DIR} || exit
eval "$(pixi shell-hook --environment prscs-construction)" # we installed the prs-cs environment with pixi
pixi info

echo "Running PRS-CS" # from https://github.com/getian107/PRScs
echo "Set parameters for efficiency [MKL_NUM_THREADS; NUMEXPR_NUM_THREADS; OMP_NUM_THREADS] ..."
export MKL_NUM_THREADS=${SLURM_NPROCS}
export NUMEXPR_NUM_THREADS=${SLURM_NPROCS}
export OMP_NUM_THREADS=${SLURM_NPROCS}
echo ""

# - Optional -
# We wanted to evaluate different parameters for phi here in the 
# beginning and included this portion to share how to do that ...

echo "Performing SLURM Array ID to Parameter Lookup ..."
PARAMETERS=$( \
    awk -v ArrayTaskID="${SLURM_ARRAY_TASK_ID}" \
    '$1 == ArrayTaskID { print $2" "$3" "$4 }' \
    "${ARRAY_CONFIG_FILE}"
)
echo "Pulled Parameter Line : { ${PARAMETERS} }"

CHROMOSOME=$( echo "$PARAMETERS" | tr ' ' '\t' | cut -f1 )
PHI_VALUE=$( echo "$PARAMETERS" | tr ' ' '\t' | cut -f2 )
POPULATION=$( echo "$PARAMETERS" | tr ' ' '\t' | cut -f3 )

echo ""
echo "!!! Running PRS-CS with [ ${POPULATION} ] on [ chr${CHROMOSOME} ] for [ PHI=${PHI_VALUE} ] !!!"
echo ""

# GWAS summary statistics formatted for PRS-CS (columns: SNP A1 A2 BETA SE)
# along with the sample size for your GWAS summary statistics (total cases and controls)
INPUT_GWAS_SUMSTATS_FILE=${PROJECT_DIR}/outputs/path/to/metaanalysis.prscs-formatted.txt
INPUT_GWAS_SAMPLE_SIZE=######

INPUT_BIM_PREFIX=${INPUTS_DIR}/cohort-name.chr${CHROMOSOME} # such as cohort1.merged.cohort2
OUTPUT_DIR=${OUTPUTS_DIR}/prscs_cohort-name
mkdir -p ${OUTPUT_DIR}

# - Optional -
# You can also add --write_psi and --write_pst to save the posterior effect
# size estimates if you have a quantitative trait and want to calculate
# individual accuracy score ...

if [ "${PHI_VALUE}" == "auto" ]; then
    echo "Running PRS-CS with automatic phi estimation ..."

    echo "+ Running PRS-CS with COHORT Files ..."
    # --ref_dir : Contains references for 1KG panels downloaded from PRScs team

    pixi run python ${PROJECT_DIR}/PRScs/PRScs.py \
        --ref_dir="${PROJECT_DIR}"/data/ldblk_"${POPULATION}" \
        --bim_prefix="${INPUT_BIM_PREFIX}" \
        --sst_file="${INPUT_GWAS_SUMSTATS_FILE}" \
        --n_gwas="${INPUT_GWAS_SAMPLE_SIZE}" \
        --chrom="${CHROMOSOME}" \
        --out_dir="${OUTPUT_DIR}"/prscs_cohort-name_"${POPULATION}" \
        --seed=9500
else
    echo "Running PRS-CS with fixed phi = ${PHI_VALUE} ..."

    echo "+ Running PRS-CS with COHORT Files ..."
    # --ref_dir : Contains references for 1KG panels downloaded from PRScs team

    pixi run python ${PROJECT_DIR}/PRScs/PRScs.py \
        --ref_dir="${PROJECT_DIR}"/data/ldblk_"${POPULATION}" \
        --bim_prefix="${INPUT_BIM_PREFIX}" \
        --sst_file="${INPUT_GWAS_SUMSTATS_FILE}" \
        --n_gwas="${INPUT_GWAS_SAMPLE_SIZE}" \
        --chrom="${CHROMOSOME}" \
        --phi="${PHI_VALUE}" \
        --out_dir="${OUTPUT_DIR}"/prscs_cohort-name_"${POPULATION}" \
        --seed=9500
fi
echo ""

echo "Ending Time: $( date )"
echo "========================================================================="

echo ""
echo "Done with script !"
echo ""
