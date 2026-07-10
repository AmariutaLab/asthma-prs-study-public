#!/bin/bash
# =============================================================================
# run_demo.sh — self-contained smoke test for the FUSION TWAS step (Stage 2)
# =============================================================================
#
# Runs FUSION.assoc_test.R on the committed demo fixture (1 GTEx tissue,
# chromosome 1, 12 genes) so anyone can confirm the TWAS association step runs
# without the full 12 GB GTEx weight panel or the controlled-access cohorts.
#
# This exercises the tissue track only. The 17-cell-type (OneK1K) track is
# PENDING — its eQTL weight panels are not publicly redistributable here (see
# README, "Cell-type (OneK1K) track — pending").
#
# ---------------------------------------------------------------------------
# DEMO DATA (see sample_data/README.md for full provenance)
#   sample_data/weights/META_Whole_Blood/*.wgt.RDat  REAL GTEx v8 EUR FUSION
#       weights, Whole_Blood, 12 genes on chr1 (11.5 KB each).
#   sample_data/LDREF/1000G.EUR.1.{bed,bim,fam}       REAL public 1000G EUR LD
#       reference, chr1 0-3.5 Mb window (978 SNPs).
#   sample_data/demo_sumstats_chr1.txt                Asthma meta-GWAS summary
#       stats (SNP Z A2 A1) from 01_meta-analysis, sliced to the demo SNPs.
#   sample_data/twas_whole_blood.demo.pos             FUSION .pos over the 12 genes.
#
# ---------------------------------------------------------------------------
# REQUIREMENTS
#   R with packages: plink2R (github: gabraham/plink2R), optparse, methods
#   FUSION TWAS toolkit (gusevlab/fusion_twas):
#       git clone https://github.com/gusevlab/fusion_twas
#   Point FUSION_DIR at that clone (or edit the default below).
#   Override R with the RSCRIPT env var if the default lacks plink2R.
# =============================================================================

set -e -o pipefail

RSCRIPT="${RSCRIPT:-Rscript}"
FUSION_DIR="${FUSION_DIR:-$HOME/fusion_twas}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SD="${SCRIPT_DIR}/sample_data"
OUT="${SCRIPT_DIR}/demo_run"
rm -rf "${OUT}"; mkdir -p "${OUT}"

echo "=========================================================================="
echo " FUSION TWAS — demo smoke test (tissue track: Whole_Blood, chr1, 12 genes)"
echo "   Rscript    : $(command -v "${RSCRIPT}" || echo NOT FOUND)"
echo "   FUSION_DIR : ${FUSION_DIR}"
echo "=========================================================================="

if [ ! -f "${FUSION_DIR}/FUSION.assoc_test.R" ]; then
    echo "ERROR: FUSION.assoc_test.R not found under FUSION_DIR=${FUSION_DIR}"
    echo "       git clone https://github.com/gusevlab/fusion_twas and set FUSION_DIR."
    exit 1
fi

cd "${SD}"
"${RSCRIPT}" "${FUSION_DIR}/FUSION.assoc_test.R" \
    --sumstats   demo_sumstats_chr1.txt \
    --weights    twas_whole_blood.demo.pos \
    --weights_dir weights \
    --ref_ld_chr LDREF/1000G.EUR. \
    --chr 1 \
    --out "${OUT}/twas_whole_blood_chr1.dat"

echo; echo ">>> TWAS result (${OUT}/twas_whole_blood_chr1.dat):"
cut -f3,6,17,19,20 "${OUT}/twas_whole_blood_chr1.dat" | column -t

echo; echo ">>> Comparing TWAS.Z against the committed reference run ..."
if command -v python3 >/dev/null; then
python3 - "$OUT/twas_whole_blood_chr1.dat" "$SD/expected_outputs/twas_whole_blood_chr1.dat" <<'PY'
import sys, csv
def z(p):
    r={}
    for row in csv.DictReader(open(p), delimiter='\t'):
        try: r[row['ID']]=float(row['TWAS.Z'])
        except: r[row['ID']]=None
    return r
got, exp = z(sys.argv[1]), z(sys.argv[2])
bad=[g for g in exp if got.get(g) is not None and exp[g] is not None and abs(got[g]-exp[g])>1e-2]
print(f"    {len(exp)} genes; mismatches (>1e-2): {len(bad)}")
print("    OK — matches reference run." if not bad else f"    DIFFERS: {bad}")
PY
fi
echo; echo "Demo complete. Output in ${OUT}/ ."
