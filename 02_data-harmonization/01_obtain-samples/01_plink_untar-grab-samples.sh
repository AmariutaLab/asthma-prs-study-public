#!/bin/bash

# This script can be used to obtain the samples present within either the GACRS
# or CAMP cohorts (i.e., after downloading the *.tar.gz files from dbGaP). We
# are using the following cohort information:

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

# For the samples we're interested in, we also ensure that IID and FID 
# information are properly put into their PLINK columns; that we have indiviudal
# sex information in the files; and that we have parental information present

################################################################################
# Note: For these variables in the script, a lot of the paths and names here
#       will be different based on your environment setup! We've removed a lot
#       of the lab-specific path information, but the code logic is the same ...
################################################################################

PROJECT_DIR=/path/to/user/project
JOB_SCRATCH=/scratch/$USER/job_$SLURM_JOB_ID # directory for temporary files on EXPANSE

COHORT_NAME=gacrs # or camp, gtex, onek1k
COHORT_NAME_UC=$( echo $COHORT_NAME | tr '[:lower:]' '[:upper:]' ) # for print statements

# Provide the full path to the downloaded dbGaP *.tar.gz archive and the 
# subdirectory name for the un-tarred *.vcf files that are generated
INPUT_TAR=${PROJECT_DIR}/data/cohort-data/name-of-tar-archive-to-process # .tar.gz
UNTAR_SUBDIR=${PROJECT_DIR}/path/to/user/vcf/output

# Prefix for PLINK output files (per chromosome), as well as the directory for
# the update files for modifying IIDs, FIDs, parent information, and sex information
PLINK_OUTPUT_PREFIX=GACRS_phs001365_TOPMed_WGS_freeze.10b.chr # example for gacrs
REFERENCE_DIR=${PROJECT_DIR}/path/to/cohort-data/plink-update-information

OUTPUTS_DIR_UNTAR=${PROJECT_DIR}/data/${UNTAR_SUBDIR} # where to store outputs
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p $OUTPUTS_DIR $OUTPUTS_DIR_UNTAR

################################################################################
echo "Untarring RAW $COHORT_NAME_UC *.tar.gz archive ..."
tar -zxvf "${INPUT_TAR}" --directory "${JOB_SCRATCH}"
echo ""

################################################################################
echo "Transferring Files from DECOMPRESSED RAW $COHORT_NAME_UC *.tar.gz archive to LUSTRE Storage ..."

# rsync options:
# -a : archive (recursive, preserves symlinks, timestamps, permissions)
# -z : compress data in transit
# -P : show progress, keep partial transfers
# -v / --stats : verbose output for debugging

rsync -azvP "${JOB_SCRATCH}"/vcf/ $OUTPUTS_DIR_UNTAR
echo ""

################################################################################
echo "Processing RAW $COHORT_NAME_UC plink files and Performing Filtering on SAMPLES OF INTEREST ..."
echo ""
echo "-------------------------------------------------------------------------"

for CHROMOSOME in $( seq 1 22 ); do
    echo "+ Running plink --update-ids for RAW ${COHORT_NAME_UC} plink file for chr${CHROMOSOME} ..."
    VCF_INPUT=$( find "$OUTPUTS_DIR_UNTAR"/*chr"${CHROMOSOME}".*.vcf.gz )
    ${PROJECT_DIR}/./plink \
        --vcf "${VCF_INPUT}" \
        --update-ids ${REFERENCE_DIR}/${COHORT_NAME}_update-plink-ids.txt \
        --make-bed \
        --out "${JOB_SCRATCH}/tmp_chr${CHROMOSOME}.update-ids"

    echo "+ Running plink --update-parents --update-sex --make-pheno for RAW ${COHORT_NAME_UC} plink file for chr${CHROMOSOME} ..."
    ${PROJECT_DIR}/./plink \
        --bfile "${JOB_SCRATCH}/tmp_chr${CHROMOSOME}.update-ids" \
        --update-sex ${REFERENCE_DIR}/${COHORT_NAME}_update-plink-sex.txt \
        --update-parents ${REFERENCE_DIR}/${COHORT_NAME}_update-plink-parents.txt \
        --must-have-sex \
        --list-duplicate-vars suppress-first \
        --make-pheno ${REFERENCE_DIR}/${COHORT_NAME}_update-plink-phenotype.txt '*' \
        --make-bed \
        --out "${OUTPUTS_DIR}/${PLINK_OUTPUT_PREFIX}${CHROMOSOME}.updated"
    rm "${JOB_SCRATCH}"/tmp_chr"${CHROMOSOME}".update-ids*

    echo "+ Running plink --keep to filter RAW ${COHORT_NAME_UC} plink file for chr${CHROMOSOME} for SAMPLES OF INTEREST ..."
    ${PROJECT_DIR}/./plink \
        --bfile "${OUTPUTS_DIR}/${PLINK_OUTPUT_PREFIX}${CHROMOSOME}.updated" \
        --keep ${REFERENCE_DIR}/${COHORT_NAME}_cases-and-controls.txt \
        --exclude "${OUTPUTS_DIR}/${PLINK_OUTPUT_PREFIX}${CHROMOSOME}.updated.dupvar" \
        --freq --nonfounders \
        --make-bed \
        --out "${OUTPUTS_DIR}/${PLINK_OUTPUT_PREFIX}${CHROMOSOME}.updated.cases-unrel-controls-only"
done

echo "-------------------------------------------------------------------------"
echo ""
echo "Done Processing RAW ${COHORT_NAME_UC} plink files and Performing Filtering on SAMPLES OF INTEREST !"
echo "Done with script !"
echo ""
