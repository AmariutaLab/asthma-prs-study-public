#!/bin/bash

# In this job, we will be doing some quality control steps with bcftools to
# handle multiallelic sites, as well as ensure that we have alligned variants to
# our reference FASTA file (i.e., from UCSC Genome Browser for our genome build)

################################################################################
#SBATCH --job-name="bcftools_clean-chromosomes"
#SBATCH --chdir="/path/to/logs/bcftools_clean-chromosomes"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=######
#SBATCH --partition=######
#SBATCH --nodes=1
#SBATCH --array=1-22
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

# Path to the indexed Reference Genome FASTA (must be indexed) along with the
# directory for the update files we used in the previous step; 
# In our script for gacrs, we used the GRCh38 reference file!
FASTA_REF_FILE=${PROJECT_DIR}/data/GRCh38_full_analysis_set_plus_decoy_hla.fa
REFERENCE_DIR=${PROJECT_DIR}/path/to/cohort-data/plink-update-information
SAMPLE_INFO=${REFERENCE_DIR}/${COHORT_NAME}_cases-and-controls.txt

# Prefix for PLINK output files (per chromosome) from the previous step, 
# but without the chromosome number, as well as the suffix to append after the
# chromosome number
PLINK_INPUT_PREFIX=GACRS_phs001365_TOPMed_WGS_freeze.10b.chr # example for gacrs
PLINK_INPUT_SUFFIX=.updated.cases-unrel-controls-only

INPUT_DIR=${PROJECT_DIR}/path/to/user/output/from/script/01
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p $OUTPUTS_DIR

################################################################################
CHROMOSOME=$SLURM_ARRAY_TASK_ID
PLINK_FILE=${INPUT_DIR}/${PLINK_INPUT_PREFIX}${CHROMOSOME}${PLINK_INPUT_SUFFIX}
echo "Processing File [ $PLINK_FILE ] ..."

################################################################################
echo "Generating List of Samples to Reorder To for ${COHORT_NAME_UC} for chr${CHROMOSOME} ..."
cut -f1,2 ${SAMPLE_INFO} | tr ' ' '_' | sort -V > "${JOB_SCRATCH}"/tmp_samples.txt
echo "Done Generating List of Samples to Reorder To for ${COHORT_NAME_UC} for chr${CHROMOSOME} !"
echo ""

echo "Processing ${COHORT_NAME_UC} File for chr${CHROMOSOME} ..."

PLINK_BASENAME=$( basename "${PLINK_FILE}" )

echo "Recoding [ ${PLINK_BASENAME} ] as VCF File instead of PLINK ..."
${PROJECT_DIR}/./plink2 \
    --bfile "${PLINK_FILE}" \
    --export vcf --output-chr chrM \
    --out "${JOB_SCRATCH}"/"${PLINK_BASENAME}".unlifted

echo "Extracting Samples from [ ${PLINK_BASENAME} ] VCF File ..."
${PROJECT_DIR}/./bcftools view \
    --samples-file "${JOB_SCRATCH}"/tmp_samples.txt \
    "${JOB_SCRATCH}"/"${PLINK_BASENAME}".unlifted.vcf \
    > "${JOB_SCRATCH}"/"${PLINK_BASENAME}".vcf

################################################################################
echo "Repairing multiallelic sites with [ bcftools norm ] ..."
${PROJECT_DIR}/./bcftools norm \
    --multiallelics -any "${JOB_SCRATCH}"/"${PLINK_BASENAME}".vcf \
    --old-rec-tag OLD_MULTIALLELIC_VAR \
    --output "${JOB_SCRATCH}"/"${PLINK_BASENAME}".bcf
${PROJECT_DIR}/./bcftools sort "${JOB_SCRATCH}"/"${PLINK_BASENAME}".bcf \
    --output "${JOB_SCRATCH}"/"${PLINK_BASENAME}".norm.bcf \
    --write-index

################################################################################
echo "Repairing alleles with [ bcftools +fixref ] ..."

# --discard      : discard sites which could not be resolved (-d)
# --fasta-ref    : reference sequence (-f)
# --mode flip    : swap or flip REF/ALT columns and GTs for non-ambiguous SNPs and ignore the rest (-m)

${PROJECT_DIR}/./bcftools +fixref \
    "${JOB_SCRATCH}"/"${PLINK_BASENAME}".norm.bcf \
    --output "${JOB_SCRATCH}"/"${PLINK_BASENAME}".bcf \
    -- \
    --discard \
    --fasta-ref ${FASTA_REF_FILE} \
    --mode flip
${PROJECT_DIR}/./bcftools sort "${JOB_SCRATCH}"/"${PLINK_BASENAME}".bcf \
    --output "${OUTPUTS_DIR}"/"${PLINK_BASENAME}".sorted.bcf \
    --write-index

# Cleanup
rm "${JOB_SCRATCH}"/tmp_samples.txt
rm "${JOB_SCRATCH}"/"${PLINK_BASENAME}".unlifted.vcf*
rm "${JOB_SCRATCH}"/"${PLINK_BASENAME}".vcf*
rm "${JOB_SCRATCH}"/"${PLINK_BASENAME}".norm.bcf*
rm "${JOB_SCRATCH}"/"${PLINK_BASENAME}".bcf*

echo "Done Processing ${COHORT_NAME_UC} File for chr${CHROMOSOME} !"
echo ""
echo "-------------------------------------------------------------------------"
echo ""
echo "Done with script !"
echo ""
