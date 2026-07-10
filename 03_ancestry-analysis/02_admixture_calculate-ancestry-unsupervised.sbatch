#!/bin/bash

################################################################################
#SBATCH --job-name="admixture_calculate-ancestry-unsupervised"
#SBATCH --chdir="/path/to/logs/admixture_calculate-ancestry-unsupervised"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --nodes=1
#SBATCH --array=1-5
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
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
echo ""

# End script automatically with errors (-e) and return error code (-o pipefail)
set -e -o pipefail
################################################################################

PROJECT_DIR=/path/to/user/project

COHORT_NAME=gacrs # or camp, gtex, onek1k
INPUT_PLINK_BED=${PROJECT_DIR}/data/cohort-data/name-of-plink-file-to-process.bed

OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p ${OUTPUTS_DIR}

################################################################################
BASENAME=$( basename "${INPUT_PLINK_BED}" )
K_VALUE=${SLURM_ARRAY_TASK_ID}

echo "Evaluating COHORT         : ${COHORT_NAME}"
echo "Running ADMIXTURE on FILE : ${BASENAME}"
echo "Running ADMIXTURE with K  : K=${K_VALUE}"
echo ""

cd ${OUTPUTS_DIR} # one quirk of admixture is it'll save outputs in the cwd
${PROJECT_DIR}/./admixture --cv "${INPUT_PLINK_BED}" "${K_VALUE}" -j"${SLURM_NTASKS_PER_NODE}"
echo ""

################################################################################
echo ""
echo "Done with script !"
echo ""
