#!/bin/bash
################################################################################
# Per-chromosome PTRS scoring driver -- tissue (39 GTEx bulk tissues),
# TWAS P+T (p-value-thresholding) variant. Identical to ptrs_score_v4.sbatch
# except it pulls LD-clumped, per-pval gene lists from
# TWAS/gene_pt_new/twas_gene_lists/ instead of the FOCUS fine-mapped lists,
# and writes to its own outputs_39_pt_* directories. For the CT version, use
# the v5 sibling scripts (ptrs_score_v5_twas_pt.sbatch +
# ptrs_score_groupC_v5_twas_pt.R).
#
# For each tissue in TISSUES, this job invokes `ptrs_score_groupC_v4_twas_pt.R`
# on the assigned chromosome to compute per-gene imputed-expression scores via
# PLINK --score and accumulate them into per-(chr, tissue)
# `all_data{,_keep,_sign}.rds` files that the concat RMDs in this folder
# consume.
#
# Configuration matrix (set COHORT below):
#
#   COHORT       OUTPUT_DIR                       PLINK basename
#   ----------   ------------------------------   -----------------------------
#   gacrs_1kg    outputs_39_pt_gacrs_1kg_new      GACRS_1kg/GACRS_1kg_chr
#   camp_1k1k    outputs_39_pt_camp_1k1k_new      CAMP_onek1k/CAMP_onek1k_chr
################################################################################
#SBATCH --job-name="ptrs_score_tissue_pt"
#SBATCH --chdir="/expanse/lustre/projects/ddp412/n5huang/gene_model/score_function/"
#SBATCH --output="%u.%x.%j.%a.stdout"
#SBATCH --error="%u.%x.%j.%a.stderr"
#SBATCH --account=csd832
#SBATCH --partition=shared
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --time=48:00:00
#SBATCH --export=ALL
################################################################################

set -euo pipefail

# ---- USER CONFIG --------------------------------------------------------------
COHORT="${COHORT:-gacrs_1kg}"           # "gacrs_1kg" or "camp_1k1k"
PVAL="${PVAL:-1}"                       # FDR/p-value tag used in output paths

PROJECT_ROOT="/expanse/lustre/projects/ddp412/n5huang"
SCRIPT_DIR="${PROJECT_ROOT}/gene_model/score_function"
# -------------------------------------------------------------------------------

echo "-------------------------------------------------------------------------"
echo "Initializing environment with PIXI"
echo "-------------------------------------------------------------------------"
cd "${SCRIPT_DIR}"
eval $(pixi shell-hook); pixi info
R --version

CHROMOSOME="${SLURM_ARRAY_TASK_ID}"

# 39 GTEx bulk tissues (tissue path table lives in TWAS/rebuilt/, 6 columns)
TWAS_PATH_TABLE="${PROJECT_ROOT}/TWAS/rebuilt/twas_path_table_39_new.csv"
GENE_LIST_TEMPLATE="${PROJECT_ROOT}/TWAS/gene_pt_new/twas_gene_lists/chr{chr}_{tissue}_gene_list_p{pval}.txt"
TISSUES=(
    Thyroid Pancreas Spleen Testis Uterus Adrenal_Gland Esophagus_Muscularis
    Esophagus_Mucosa Cells_Cultured_fibroblasts Skin_Sun_Exposed_Lower_leg
    Artery_Tibial Nerve_Tibial Heart_Atrial_Appendage Brain_Spinal_cord_cervical_c-1
    Liver Cells_EBV-transformed_lymphocytes Lung Prostate Muscle_Skeletal Ovary
    Brain_Cortex Brain_BasalGanglia Brain_Cerebellum Whole_Blood Adipose_Subcutaneous
    Adipose_Visceral_Omentum Artery_Aorta Colon_Sigmoid Colon_Transverse Brain_Limbic
    Artery_Coronary Brain_Substantia_nigra Breast_Mammary_Tissue Heart_Left_Ventricle
    Minor_Salivary_Gland Pituitary Skin_Not_Sun_Exposed_Suprapubic Stomach Vagina
)

case "${COHORT}" in
    gacrs_1kg)
        PLINK_BASENAME="${PROJECT_ROOT}/GACRS_1kg/GACRS_1kg_chr"
        OUTPUT_DIR_BASE="${SCRIPT_DIR}/outputs_39_pt_gacrs_1kg_new"
        ;;
    camp_1k1k)
        PLINK_BASENAME="${PROJECT_ROOT}/CAMP_onek1k/CAMP_onek1k_chr"
        OUTPUT_DIR_BASE="${SCRIPT_DIR}/outputs_39_pt_camp_1k1k_new"
        ;;
    *)
        echo "ERROR: unknown COHORT='${COHORT}' (expected 'gacrs_1kg' or 'camp_1k1k')" >&2
        exit 1
        ;;
esac

echo "FEATURE_SET=tissue  COHORT=${COHORT}  PVAL=${PVAL}  CHR=${CHROMOSOME}"
echo "TWAS_PATH_TABLE=${TWAS_PATH_TABLE}"
echo "PLINK_BASENAME=${PLINK_BASENAME}"
echo "GENE_LIST_TEMPLATE=${GENE_LIST_TEMPLATE}"
echo "OUTPUT_DIR_BASE=${OUTPUT_DIR_BASE}"

echo "-------------------------------------------------------------------------"
echo "Running ptrs_score_groupC_v4_twas_pt.R per tissue"
echo "-------------------------------------------------------------------------"

for TISSUE in "${TISSUES[@]}"; do
    echo "-------------------------------------------------------------------------"
    echo ">>>>> Evaluating ${TISSUE}"

    TISSUE_OUTPUT_DIR="${OUTPUT_DIR_BASE}/${PVAL}/${TISSUE}"
    mkdir -p "${TISSUE_OUTPUT_DIR}"

    # Tissue TWAS path table has 6 columns:
    #   tissue, gtex_twas, transcripts, transcript_keep_file, weights, PANEL
    # transcript_keep_file is a real subset filter (e.g. heritable-genes list).
    line=$(grep "^\"${TISSUE}\"," "${TWAS_PATH_TABLE}" || true)
    if [[ -z "${line}" ]]; then
        echo "Error: no entry for tissue '${TISSUE}' in ${TWAS_PATH_TABLE}" >&2
        exit 1
    fi
    line_clean=$(echo "${line}" | tr -d '"')
    IFS=',' read -r tissue_name gtex_file transcript_file transcript_filter_file weight_file_bname panel <<< "${line_clean}"

    echo "---- PARSED VARIABLES ----"
    echo "tissue_name=${tissue_name}"
    echo "gtex_file=${gtex_file}"
    echo "transcripts=${transcript_file}"
    echo "filter_file=${transcript_filter_file}"
    echo "weights=${weight_file_bname}"
    echo "PANEL=${panel}"
    echo "--------------------------"

    Rscript "${SCRIPT_DIR}/ptrs_score_groupC_v4_twas_pt.R" \
        "${CHROMOSOME}" \
        "${transcript_file}" \
        "${gtex_file}" \
        "${transcript_filter_file}" \
        "${weight_file_bname}" \
        "${TISSUE_OUTPUT_DIR}" \
        "${TISSUE}" \
        "${PVAL}" \
        "${PLINK_BASENAME}" \
        "${GENE_LIST_TEMPLATE}"
done

echo "-------------------------------------------------------------------------"
echo "Done!"
echo "========================================================================="
