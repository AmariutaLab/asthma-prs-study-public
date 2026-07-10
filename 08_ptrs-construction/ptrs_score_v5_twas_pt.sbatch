#!/bin/bash
################################################################################
# Per-chromosome PTRS scoring driver -- CT (17 OneK1K cell types),
# TWAS P+T (p-value-thresholding) variant. Identical to ptrs_score_v5.sbatch
# except it pulls LD-clumped, per-pval gene lists from
# TWAS/gene_pt_ct_new/twas_gene_lists/ instead of the FOCUS fine-mapped lists,
# and writes to its own outputs_17CT_pt_* directories. For the tissue version,
# use the v4 sibling scripts (ptrs_score_v4_twas_pt.sbatch +
# ptrs_score_groupC_v4_twas_pt.R).
#
# For each cell type in TISSUES, this job invokes
# `ptrs_score_groupC_v5_twas_pt.R` on the assigned chromosome to compute
# per-gene imputed-expression scores via PLINK --score and accumulate them
# into per-(chr, cell-type) `all_data{,_keep,_sign}.rds` files that the concat
# RMDs in this folder consume.
#
# Configuration matrix (set COHORT below):
#
#   COHORT       OUTPUT_DIR                     PLINK basename
#   ----------   ----------------------------   -----------------------------
#   gacrs_1kg    outputs_17CT_pt_gacrs_1kg_new  GACRS_1kg/GACRS_1kg_chr
#   camp_gtex    outputs_17CT_pt_camp_gtex_new  CAMP_GTEX/CAMP_GTEX_chr
################################################################################
#SBATCH --job-name="ptrs_score_ct_pt"
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
COHORT="${COHORT:-camp_gtex}"           # "gacrs_1kg" or "camp_gtex"
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

# 17 OneK1K cell types (CT-only path tables live in TWAS/rebuilt_ct/, 5 columns)
TWAS_PATH_TABLE="${PROJECT_ROOT}/TWAS/rebuilt_ct/twas_path_table_ct_new.csv"
GENE_LIST_TEMPLATE="${PROJECT_ROOT}/TWAS/gene_pt_ct_new/twas_gene_lists/chr{chr}_{tissue}_gene_list_p{pval}.txt"
TISSUES=(
    b_intermediate b_memory b_naive cd14_mono cd16_mono cd4_ctl cd4_naive cd4_tcm
    cd4_tem cd8_naive cd8_tcm cd8_tem gdt mait nk nk_cd56bright treg
)

case "${COHORT}" in
    gacrs_1kg)
        PLINK_BASENAME="${PROJECT_ROOT}/GACRS_1kg/GACRS_1kg_chr"
        OUTPUT_DIR_BASE="${SCRIPT_DIR}/outputs_17CT_pt_gacrs_1kg_new"
        ;;
    camp_gtex)
        PLINK_BASENAME="${PROJECT_ROOT}/CAMP_GTEX/CAMP_GTEX_chr"
        OUTPUT_DIR_BASE="${SCRIPT_DIR}/outputs_17CT_pt_camp_gtex_new"
        ;;
    *)
        echo "ERROR: unknown COHORT='${COHORT}' (expected 'gacrs_1kg' or 'camp_gtex')" >&2
        exit 1
        ;;
esac

echo "FEATURE_SET=ct  COHORT=${COHORT}  PVAL=${PVAL}  CHR=${CHROMOSOME}"
echo "TWAS_PATH_TABLE=${TWAS_PATH_TABLE}"
echo "PLINK_BASENAME=${PLINK_BASENAME}"
echo "GENE_LIST_TEMPLATE=${GENE_LIST_TEMPLATE}"
echo "OUTPUT_DIR_BASE=${OUTPUT_DIR_BASE}"

echo "-------------------------------------------------------------------------"
echo "Running ptrs_score_groupC_v5_twas_pt.R per cell type"
echo "-------------------------------------------------------------------------"

for TISSUE in "${TISSUES[@]}"; do
    echo "-------------------------------------------------------------------------"
    echo ">>>>> Evaluating ${TISSUE}"

    TISSUE_OUTPUT_DIR="${OUTPUT_DIR_BASE}/${PVAL}/${TISSUE}"
    mkdir -p "${TISSUE_OUTPUT_DIR}"

    # CT TWAS path table has 5 columns:
    #   tissue, gtex_twas, transcripts, weights, PANEL
    # No separate filter column -- pass transcripts as the (no-op) filter.
    line=$(grep "^\"${TISSUE}\"," "${TWAS_PATH_TABLE}" || true)
    if [[ -z "${line}" ]]; then
        echo "Error: no entry for tissue '${TISSUE}' in ${TWAS_PATH_TABLE}" >&2
        exit 1
    fi
    line_clean=$(echo "${line}" | tr -d '"')
    IFS=',' read -r tissue_name gtex_file transcript_file weight_file_bname panel <<< "${line_clean}"
    transcript_filter_file="${transcript_file}"

    echo "---- PARSED VARIABLES ----"
    echo "tissue_name=${tissue_name}"
    echo "gtex_file=${gtex_file}"
    echo "transcripts=${transcript_file}"
    echo "filter_file=${transcript_filter_file}"
    echo "weights=${weight_file_bname}"
    echo "PANEL=${panel}"
    echo "--------------------------"

    Rscript "${SCRIPT_DIR}/ptrs_score_groupC_v5_twas_pt.R" \
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
