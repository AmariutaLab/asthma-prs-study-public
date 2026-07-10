#!/bin/bash
# =============================================================================
# run_demo.sh — self-contained smoke test for the ADMIXTURE ancestry pipeline
# =============================================================================
#
# This is a small, non-SLURM demonstration of the two analysis steps documented
# in this directory's README, wired up to run end-to-end on the committed demo
# data in sample_data/.  It exists so anyone who clones the repo can confirm the
# pipeline logic runs without access to the controlled-access study cohorts.
#
#   Step 01  ->  plink : filter to HapMap3 variants + variant QC + LD pruning
#   Step 02  ->  admixture --cv : unsupervised ancestry for K = 1..3
#
# The canonical, cluster-oriented versions of these steps are:
#   01_plink_filter-hm3-variants.sh
#   02_admixture_calculate-ancestry-unsupervised.sbatch
# run_demo.sh mirrors their logic but points at the demo inputs and drops the
# SLURM / lab-path scaffolding.
#
# ---------------------------------------------------------------------------
# DEMO DATA  (see sample_data/README.md for full provenance)
#   sample_data/demo_1kg.chr22.{bed,bim,fam}  REAL public 1000 Genomes Project
#       Phase 3 EUR genotypes, chr22 only, downsampled to 100 samples / 6,000
#       SNPs.  NOT individual-level study data.
#   sample_data/w_hm3_chr22_subset.snplist    REAL HapMap3 SNP list (Alkes Group,
#       Zenodo 7773502), subset to the chr22 rsIDs present in the demo.
#
# NOTE ON INTERPRETATION: the demo genotypes are single-ancestry (1000G EUR),
# so this is a *smoke test of the pipeline*, not a biological result — the
# cross-validation error will simply favour K=1.  With a genuinely multi-ancestry
# input (e.g. study cohort merged with the full 1000 Genomes reference) the same
# commands surface real population structure.
#
# ---------------------------------------------------------------------------
# REQUIREMENTS
#   plink      (PLINK 1.9)      https://www.cog-genomics.org/plink/
#   admixture  (ADMIXTURE 1.3)  https://dalexander.github.io/admixture/
# Both are on bioconda:  conda create -n asthma-demo -c bioconda plink admixture
# Override the binaries with the PLINK / ADMIXTURE env vars if not on PATH.
# =============================================================================

set -e -o pipefail

PLINK="${PLINK:-plink}"
ADMIXTURE="${ADMIXTURE:-admixture}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SAMPLE_DIR="${SCRIPT_DIR}/sample_data"
OUTPUTS_DIR="${SCRIPT_DIR}/demo_run"
rm -rf "${OUTPUTS_DIR}"; mkdir -p "${OUTPUTS_DIR}"

INPUT_PLINK="${SAMPLE_DIR}/demo_1kg.chr22"
HAPMAP3_SNPLIST="${SAMPLE_DIR}/w_hm3_chr22_subset.snplist"
BASENAME="demo_1kg.chr22"

echo "=========================================================================="
echo " ADMIXTURE ancestry pipeline — demo smoke test"
echo "   plink     : $(command -v "${PLINK}" || echo "NOT FOUND")"
echo "   admixture : $(command -v "${ADMIXTURE}" || echo "NOT FOUND")"
echo "   input     : ${INPUT_PLINK}.{bed,bim,fam}"
echo "   outputs   : ${OUTPUTS_DIR}"
echo "=========================================================================="

# --- Step 01 : filter to HapMap3 variants, QC, LD prune ---------------------
echo; echo ">>> Step 01 — filter to HapMap3 variants + QC + LD pruning"

# strip the header off the HapMap3 snplist (col 1 = rsID)
awk 'NR > 1 {print $1}' "${HAPMAP3_SNPLIST}" > "${OUTPUTS_DIR}/tmp_hm3_rsids.txt"

"${PLINK}" --bfile "${INPUT_PLINK}" \
    --extract "${OUTPUTS_DIR}/tmp_hm3_rsids.txt" \
    --allow-no-sex \
    --make-bed --out "${OUTPUTS_DIR}/${BASENAME}"

"${PLINK}" --bfile "${OUTPUTS_DIR}/${BASENAME}" \
    --geno 0.05 --maf 0.01 \
    --make-bed --out "${OUTPUTS_DIR}/${BASENAME}.variant-qc"

"${PLINK}" --bfile "${OUTPUTS_DIR}/${BASENAME}.variant-qc" \
    --indep-pairwise 50 10 0.1 \
    --out "${OUTPUTS_DIR}/${BASENAME}.ld-qc"

"${PLINK}" --bfile "${OUTPUTS_DIR}/${BASENAME}.variant-qc" \
    --extract "${OUTPUTS_DIR}/${BASENAME}.ld-qc.prune.in" \
    --make-bed --out "${OUTPUTS_DIR}/${BASENAME}.filt-hm3"

N_SNP=$( wc -l < "${OUTPUTS_DIR}/${BASENAME}.filt-hm3.bim" )
echo "    -> ${BASENAME}.filt-hm3 : ${N_SNP} variants retained"

# --- Step 02 : unsupervised ADMIXTURE for K = 1..3 --------------------------
echo; echo ">>> Step 02 — unsupervised ADMIXTURE (K = 1..3)"
cd "${OUTPUTS_DIR}"                       # ADMIXTURE writes to the cwd
: > cv_error.txt
for K in 1 2 3; do
    echo "    running ADMIXTURE K=${K} ..."
    "${ADMIXTURE}" --cv "${BASENAME}.filt-hm3.bed" "${K}" -j2 \
        | tee "admixture.K${K}.log" \
        | grep -i "CV error" >> cv_error.txt || true
done

echo; echo ">>> Cross-validation error by K:"; cat cv_error.txt

echo; echo "=========================================================================="
echo " Demo complete. Outputs in ${OUTPUTS_DIR}"
echo "   ${BASENAME}.filt-hm3.{bed,bim,fam}   filtered/pruned fileset"
echo "   ${BASENAME}.filt-hm3.{K}.Q           ancestry proportions per sample"
echo "   ${BASENAME}.filt-hm3.{K}.P           allele frequencies per variant"
echo "   cv_error.txt                         CV error used to select optimal K"
echo " A committed reference run is in sample_data/expected_outputs/ ."
echo "=========================================================================="
