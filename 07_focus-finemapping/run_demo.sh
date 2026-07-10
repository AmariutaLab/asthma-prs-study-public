#!/bin/bash
# =============================================================================
# run_demo.sh — self-contained smoke test for ma-FOCUS fine-mapping (Stage 2)
# =============================================================================
#
# Runs `focus import` + `focus finemap` on the committed demo fixture
# (Whole_Blood, chr1, 12 genes) so anyone can confirm the fine-mapping step runs
# without the full 12 GB GTEx weight panel or the controlled-access cohorts.
#
# IMPORTANT: this uses the VENDORED, locally-modified ma-FOCUS shipped in
# ./ma-focus-modified/ (patched to import the meta-analyzed bulk-tissue FUSION
# `META_*` panels, which crash stock ma-focus). Install it before running:
#
#     conda create -n ma-focus python=3.8 -y
#     conda activate ma-focus
#     pip install ./ma-focus-modified
#     pip install mygene 'rpy2==3.5.12'
#
# rpy2 loads FUSION *.wgt.RDat via R, so R_HOME must point at an R whose
# architecture matches your Python (e.g. an x86_64 conda R for x86_64 Python):
#
#     conda install -n ma-focus r-base -y     # simplest: gives the env its own R
#     export R_HOME="$CONDA_PREFIX/lib/R"     # or any same-arch R install root
#
# Tissue track only. The 17-cell-type (OneK1K) track is PENDING (see README).
#
# DEMO DATA (see sample_data/README.md):
#   sample_data/weights/META_Whole_Blood/*.wgt.RDat   REAL GTEx v8 EUR FUSION weights (12 chr1 genes)
#   sample_data/twas_whole_blood.demo.pos             .pos over the 12 genes (WGT relative to sample_data/)
#   sample_data/sumstat/chr1_..._sumstats.mod2        asthma meta-GWAS in FOCUS format (SNP Z A2 A1 CHR BP N)
#   ../05_twas-analysis/sample_data/LDREF/1000G.EUR.1 REAL public 1000G EUR LD reference (chr1 window)
# =============================================================================

set -e -o pipefail

FOCUS="${FOCUS:-focus}"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SD="${SCRIPT_DIR}/sample_data"
LDREF="${LDREF:-${SCRIPT_DIR}/../05_twas-analysis/sample_data/LDREF/1000G.EUR.1}"
OUT="${SCRIPT_DIR}/demo_run"
rm -rf "${OUT}"; mkdir -p "${OUT}/results/Whole_Blood"

echo "=========================================================================="
echo " ma-FOCUS fine-mapping — demo smoke test (Whole_Blood, chr1, 12 genes)"
echo "   focus  : $(command -v "${FOCUS}" || echo NOT FOUND)"
echo "   R_HOME : ${R_HOME:-<unset — rpy2 will use the default R>}"
echo "   LDREF  : ${LDREF}"
echo "=========================================================================="
if ! command -v "${FOCUS}" >/dev/null; then
    echo "ERROR: 'focus' not found. Install the vendored ma-focus-modified (see header)."; exit 1
fi

cd "${SD}"          # so the .pos WGT paths (weights/...) resolve

echo; echo ">>> focus import (FUSION weights -> FOCUS .db)"
"${FOCUS}" import twas_whole_blood.demo.pos fusion \
    --name GTEx --assay rnaseq --output "${OUT}/fusion_demo"

echo; echo ">>> focus finemap (chr1)"
"${FOCUS}" finemap \
    "${SD}/sumstat/chr1_formatted_meta_analysis.sumstats.mod2" \
    "${LDREF}" \
    "${OUT}/fusion_demo.db" \
    --chr 1 --prior-prob gencode38 --locations 38:EUR --p-threshold 1 \
    --out "${OUT}/results/Whole_Blood/focus_result_demo"

RES="${OUT}/results/Whole_Blood/focus_result_demo.focus.tsv"
echo; echo ">>> credible-set genes (in_cred_set_pop1 == 1):"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="mol_name")m=i; if($i=="pips_pop1")p=i; if($i=="in_cred_set_pop1")c=i}}
            NR>1 && $c==1 && $2!="NULL.MODEL"{printf "    %-12s pip=%s\n", $4, $p}' "${RES}"

echo; echo ">>> comparing credible-set membership against the committed reference run ..."
EXP="${SD}/expected_outputs/focus_result_demo.focus.tsv"
cs(){ awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="mol_name")m=i; if($i=="in_cred_set_pop1")c=i}}
                  NR>1 && $c==1 && $2!="NULL.MODEL"{print $4}' "$1" | sort -u; }
if diff <(cs "${RES}") <(cs "${EXP}") >/dev/null; then echo "    OK — credible set matches reference."; else
    echo "    NOTE: credible set differs from reference:"; diff <(cs "${RES}") <(cs "${EXP}") || true; fi

echo; echo "Demo complete. Output in ${OUT}/ ."
