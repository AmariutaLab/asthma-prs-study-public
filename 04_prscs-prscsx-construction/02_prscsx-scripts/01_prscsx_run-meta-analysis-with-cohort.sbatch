#!/bin/bash

################################################################################
#SBATCH --job-name="prscsx_run-meta-analysis-with-cohort"
#SBATCH --chdir="/path/to/logs/prscsx_run-meta-analysis-with-cohort"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --array=1-66
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=36G
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
ARRAY_CONFIG_FILE=${SLURM_SUBMIT_DIR}/prscsx_run-parameters.txt

# We're using the PLINK file to be used as a validation target dataset
# containing the SNPs of interest within our summary statistics reference
INPUTS_DIR=${PROJECT_DIR}/data/cohort-data/name-of-plink-file-to-process
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p ${OUTPUTS_DIR}

################################################################################
echo "Initializing Environment with PIXI"
cd ${PROJECT_DIR} || exit
eval "$(pixi shell-hook --environment prscs-construction)" # we installed the prs-cs environment with pixi (same as prs-csx)
pixi info

echo "Running PRS-CSx" # from https://github.com/getian107/PRScsx
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
    '$1 == ArrayTaskID { print $2" "$3 }' \
    "${ARRAY_CONFIG_FILE}"
)
echo "Pulled Parameter Line : { ${PARAMETERS} }"

CHROMOSOME=$( echo "$PARAMETERS" | tr ' ' '\t' | cut -f1 )
PHI_VALUE=$( echo "$PARAMETERS" | tr ' ' '\t' | cut -f2 )

echo ""
echo "!!! Running PRS-CSx with the 1000 Genomes LD Block Panel on [ chr${CHROMOSOME} ] for [ PHI=${PHI_VALUE} ] !!!"
echo ""

################################################################################
INPUT_BIM_PREFIX=${INPUTS_DIR}/cohort-name.chr${CHROMOSOME} # such as cohort1.merged.cohort2
OUTPUT_DIR=${OUTPUTS_DIR}/prscsx_cohort-name
mkdir -p ${OUTPUT_DIR}

# EUR GWAS summary statistics formatted for PRS-CS (columns: SNP A1 A2 BETA SE)
# along with the sample size for your GWAS summary statistics (total cases and controls)
INPUT_GWAS_SUMSTATS_EUR=${PROJECT_DIR}/outputs/path/to/metaanalysis.prscsx-formatted.txt
INPUT_GWAS_SAMPLE_SIZE_EUR=######

# We also included ancestry-specific summary statistics to use multi-ancestry
# capabilities of PRS-CSx by leveraging GBMI summary statistics, as EUR
# was our meta-analysis of GBMI and TAGC
INPUT_GWAS_SUMSTATS_AFR=${PROJECT_DIR}/outputs/path/to/GBMI_sumstats_afr.prscsx-formatted.txt
INPUT_GWAS_SUMSTATS_AMR=${PROJECT_DIR}/outputs/path/to/GBMI_sumstats_amr.prscsx-formatted.txt
INPUT_GWAS_SUMSTATS_EAS=${PROJECT_DIR}/outputs/path/to/GBMI_sumstats_eas.prscsx-formatted.txt
INPUT_GWAS_SUMSTATS_SAS=${PROJECT_DIR}/outputs/path/to/GBMI_sumstats_sas.prscsx-formatted.txt

INPUT_GWAS_SAMPLE_SIZE_AFR=######  # total cases + controls for AFR GWAS
INPUT_GWAS_SAMPLE_SIZE_AMR=######  # total cases + controls for AMR GWAS
INPUT_GWAS_SAMPLE_SIZE_EAS=######  # total cases + controls for EAS GWAS
INPUT_GWAS_SAMPLE_SIZE_SAS=######  # total cases + controls for SAS GWAS

INPUT_SST_FILE=${INPUT_GWAS_SUMSTATS_EUR},${INPUT_GWAS_SUMSTATS_AFR},${INPUT_GWAS_SUMSTATS_AMR},${INPUT_GWAS_SUMSTATS_EAS},${INPUT_GWAS_SUMSTATS_SAS}
INPUT_N_GWAS=${INPUT_GWAS_SAMPLE_SIZE_EUR},${INPUT_GWAS_SAMPLE_SIZE_AFR},${INPUT_GWAS_SAMPLE_SIZE_AMR},${INPUT_GWAS_SAMPLE_SIZE_EAS},${INPUT_GWAS_SAMPLE_SIZE_SAS}

################################################################################
if [ "${PHI_VALUE}" == "auto" ]; then
    echo "Running PRS-CSx with automatic phi estimation ..."

    echo "+ Running PRS-CSx with COHORT Files ..."
    pixi run python ${PROJECT_DIR}/PRScsx/PRScsx.py \
        --ref_dir="${PROJECT_DIR}"/data/ldblk_1kg \
        --bim_prefix="${INPUT_BIM_PREFIX}" \
        --sst_file="${INPUT_SST_FILE}" \
        --n_gwas="${INPUT_N_GWAS}" \
        --pop=EUR,AFR,AMR,EAS,SAS \
        --chrom="${CHROMOSOME}" \
        --out_dir="${OUTPUT_DIR}" \
        --out_name=prscsx_cohort-name_1kg \
        --meta=True \
        --seed=9500
else
    echo "Running PRS-CSx with fixed phi = ${PHI_VALUE} ..."

    echo "+ Running PRS-CSx with COHORT Files ..."
    pixi run python ${PROJECT_DIR}/PRScsx/PRScsx.py \
        --ref_dir="${PROJECT_DIR}"/data/ldblk_1kg \
        --bim_prefix="${INPUT_BIM_PREFIX}" \
        --sst_file="${INPUT_SST_FILE}" \
        --n_gwas="${INPUT_N_GWAS}" \
        --pop=EUR,AFR,AMR,EAS,SAS \
        --chrom="${CHROMOSOME}" \
        --phi="${PHI_VALUE}" \
        --out_dir="${OUTPUT_DIR}" \
        --out_name=prscsx_cohort-name_1kg \
        --meta=True \
        --seed=9500
fi

echo ""
echo "Ending Time: $( date )"
echo "========================================================================="

echo ""
echo "Done with Script !"
echo ""
