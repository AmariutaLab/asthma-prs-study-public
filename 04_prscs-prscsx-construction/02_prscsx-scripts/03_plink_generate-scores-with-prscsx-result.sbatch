#!/bin/bash

################################################################################
#SBATCH --job-name="plink_generate-scores-with-prscsx-result"
#SBATCH --chdir="/path/to/logs/plink_generate-scores-with-prscsx-result"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-6
#SBATCH --mem=36G
#SBATCH --time=24:00:00
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

INPUTS_DIR=${PROJECT_DIR}/path/to/user/output/from/script/02
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p ${OUTPUTS_DIR}

################################################################################
echo "Initializing Environment with PIXI"
cd ${PROJECT_DIR} || exit
eval "$(pixi shell-hook --environment prscs-construction)"
pixi info

################################################################################
COHORT="cohort-name"
COHORT_PLINK_NAME="cohort1.merged.cohort2" # match per-chromosome PLINK filenames
POPULATION=""

if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    POPULATION="AFR"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 2 ]; then
    POPULATION="AMR"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 3 ]; then
    POPULATION="EAS"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 4 ]; then
    POPULATION="EUR"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 5 ]; then
    POPULATION="SAS"
elif [ "${SLURM_ARRAY_TASK_ID}" -eq 6 ]; then
    POPULATION="META"
else
    echo "Error: Invalid SLURM_ARRAY_TASK_ID value: ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Merging Per-Chromosome PLINK files to Create COHORT 1-COHORT 2 Combined Set ..."
TMP_MERGE_LIST="${JOB_SCRATCH}"/tmp_merge_"${COHORT}"_plink_list.txt
find "${PROJECT_DIR}"/data/cohort-data/dir-of-plink-file-to-process/*"${COHORT}"*.bim | \
    sort -V | sed 's/.bim//g' > "${TMP_MERGE_LIST}"
TMP_MERGE_PLINK="${JOB_SCRATCH}"/"${COHORT_PLINK_NAME}".allchr
${PROJECT_DIR}/./plink --merge-list "${TMP_MERGE_LIST}" --make-bed --out "${TMP_MERGE_PLINK}"
echo ""

echo "Scoring PLINK files for COHORT 1-COHORT 2 Combined Set with PRS-CSx Weights in ${POPULATION} Individuals ..."
for PHI_VALUE in 1e-02 1e-04 auto; do
    echo "---------------------------------------------------------------"
    echo "+ Calculating Score for PHI=${PHI_VALUE} ..."
    INPUT_PRSCSX_SCORE_FN=${INPUTS_DIR}/prscsx_${COHORT}_1kg_${POPULATION}_pst_eff_a1_b0.5_phi${PHI_VALUE}_allchr.txt
    OUTPUT_PLINK_SCORE_FN=${OUTPUTS_DIR}/${COHORT}/${COHORT}_allchr_prscsx_phi${PHI_VALUE}_1kg_${POPULATION}_score
    ${PROJECT_DIR}/./plink \
        --bfile "${TMP_MERGE_PLINK}" \
        --score "${INPUT_PRSCSX_SCORE_FN}" 2 4 6 sum center \
        --out "${OUTPUT_PLINK_SCORE_FN}"
    echo "Done Calculating Score for PHI=${PHI_VALUE} !"
done

echo "Ending Time: $( date )"
echo "=========================================================================="

echo ""
echo "Done with Script !"
echo ""
