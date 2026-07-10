#!/bin/bash
# =============================================================================
# run_demo.sh — self-contained smoke test for PTRS scoring (Stage 1, tissue)
# =============================================================================
#
# Runs the ACTUAL Stage-1 worker `ptrs_score_groupC_v4.R` (unmodified) on the
# committed demo fixture (Whole_Blood, chr1, 12 genes), so anyone can confirm
# the per-(chr, tissue) PTRS scoring runs without the full GTEx weight panel or
# the controlled-access cohorts.
#
# The worker imputes predicted expression per gene (plink --score), weights by
# the per-gene TWAS z, and sums across genes into three RDS files:
#   chr1_Whole_Blood_all_data.rds        raw sum
#   chr1_Whole_Blood_all_data_keep.rds   sum weighted by signed z-score
#   chr1_Whole_Blood_all_data_sign.rds   sum weighted by sign(z-score)
#
# The worker calls `${PROJECT_DIR}/./plink` via system(); this script points
# PROJECT_DIR at a demo_run/ dir containing a `plink` symlink, so the real
# worker runs verbatim. It also copies the genotype target into demo_run/ first,
# because the worker rewrites the .bim in place (strand/allele harmonization).
#
# Tissue track only. The 17-cell-type (OneK1K) track is PENDING (see README).
#
# DEMO DATA (see sample_data/README.md):
#   sample_data/weights/META_Whole_Blood/*.wgt.RDat        REAL GTEx v8 EUR FUSION weights (12 chr1 genes)
#   sample_data/TranscriptsInWhole_BloodModel.demo.txt     the 12 gene ids
#   sample_data/Marginal_alphas_...Whole_Blood.demo.txt.gz per-gene TWAS z (from the 05 demo)
#   sample_data/focus_geneList_39/chr1_Whole_Blood_gene_list_p1.txt  selected-gene list
#   genotype target: ../05_twas-analysis/sample_data/LDREF/1000G.EUR.1 (public 1000G EUR,
#       standing in for the controlled-access cohort target)
#
# REQUIREMENTS: R with data.table, glmnet, dplyr; PLINK 1.9 on PATH (or set PLINK).
# =============================================================================

set -e -o pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
PLINK="${PLINK:-plink}"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SD="${SCRIPT_DIR}/sample_data"
TARGET="${TARGET:-${SCRIPT_DIR}/../05_twas-analysis/sample_data/LDREF/1000G.EUR.1}"
OUT="${SCRIPT_DIR}/demo_run"
rm -rf "${OUT}"; mkdir -p "${OUT}/Whole_Blood"

PLINK_BIN="$(command -v "${PLINK}" || true)"
echo "=========================================================================="
echo " PTRS scoring — demo smoke test (Whole_Blood, chr1, 12 genes)"
echo "   Rscript : $(command -v "${RSCRIPT}" || echo NOT FOUND)"
echo "   plink   : ${PLINK_BIN:-NOT FOUND}"
echo "   target  : ${TARGET}.{bed,bim,fam}"
echo "=========================================================================="
if [ -z "${PLINK_BIN}" ]; then echo "ERROR: plink not found (set PLINK=/path/to/plink)."; exit 1; fi

# plink symlink so the worker's ${PROJECT_DIR}/./plink resolves
ln -sf "${PLINK_BIN}" "${OUT}/plink"
# copy the target genotypes (worker rewrites the .bim in place)
for e in bed bim fam; do cp "${TARGET}.${e}" "${OUT}/target_chr1.${e}"; done

export PROJECT_DIR="${OUT}"
"${RSCRIPT}" "${SCRIPT_DIR}/ptrs_score_groupC_v4.R" \
    1 \
    "${SD}/TranscriptsInWhole_BloodModel.demo.txt" \
    "${SD}/Marginal_alphas_NEW_TWAS_Whole_Blood.demo.txt.gz" \
    "${SD}/TranscriptsInWhole_BloodModel.demo.txt" \
    "${SD}/weights/META_Whole_Blood/Whole_Blood." \
    "${OUT}/Whole_Blood" \
    Whole_Blood \
    1 \
    "${OUT}/target_chr" \
    "${SD}/focus_geneList_39/chr{chr}_{tissue}_gene_list_p{pval}.txt" \
    > "${OUT}/worker.log" 2>&1 || { echo "worker failed; see ${OUT}/worker.log"; tail -20 "${OUT}/worker.log"; exit 1; }

echo; echo ">>> per-(chr,tissue) PTRS scores written:"
"${RSCRIPT}" -e '
for (f in c("all_data","all_data_keep","all_data_sign")) {
  d <- readRDS(file.path("'"${OUT}/Whole_Blood"'", paste0("chr1_Whole_Blood_", f, ".rds")))
  cat(sprintf("    chr1_Whole_Blood_%-13s %d samples x %d cols  SCORE[1:3] = %s\n",
      paste0(f,".rds"), nrow(d), ncol(d), paste(round(head(d$SCORE,3),3), collapse=", ")))
}' 2>/dev/null

echo; echo ">>> comparing SCORE against the committed reference run ..."
"${RSCRIPT}" -e '
ok <- TRUE
for (f in c("all_data","all_data_keep","all_data_sign")) {
  a <- readRDS(file.path("'"${OUT}/Whole_Blood"'", paste0("chr1_Whole_Blood_",f,".rds")))
  b <- readRDS(file.path("'"${SD}/expected_outputs"'", paste0("chr1_Whole_Blood_",f,".rds")))
  if (max(abs(a$SCORE - b$SCORE)) > 1e-6) { ok <- FALSE; cat("    DIFFERS:", f, "\n") }
}
cat(if (ok) "    OK — matches reference run.\n" else "    (differences above)\n")' 2>/dev/null

echo; echo "Demo complete. Outputs in ${OUT}/Whole_Blood/ ."
echo "Stage 2 (cross-chromosome concat + cohort split) is the *_unified.RMD notebooks;"
echo "they consume these per-chromosome RDS files — see README."
