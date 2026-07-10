#!/bin/bash
# =============================================================================
# run_demo.sh — self-contained smoke test for PRS-CS + PRS-CSx construction
# =============================================================================
#
# Runs PRS-CS and PRS-CSx end-to-end on the committed demo fixture
# (chromosome 22, 250 HapMap3 SNPs), then applies the posterior weights to a
# validation target with `plink --score` — so anyone can confirm the pipeline
# runs without the official ~4.5 GB LD reference panel or controlled cohorts.
#
#   Step 01  ->  PRS-CS / PRS-CSx : posterior SNP effect sizes (MCMC)
#   Step 02  ->  merge posteriors across chromosomes (trivial here: 1 chr)
#   Step 03  ->  plink --score : per-sample PRS on the target
#
# DEMO DATA (see sample_data/README.md for full provenance)
#   reference/ldblk_1kg_eur/{snpinfo_1kg_hm3, ldblk_1kg_chr22.hdf5}   REAL 1000G
#       EUR LD reference built from public chr22 genotypes (250 SNPs, 2 blocks).
#   reference/snpinfo_mult_1kg_hm3                                    PRS-CSx multi-pop snpinfo.
#   demo_sumstats_chr22.txt   FAKE, SIMULATED GWAS sumstats (SNP A1 A2 BETA SE) —
#       effect sizes are random draws over the demo SNPs, NOT real asthma
#       associations. Format matches the 01_meta-analysis output. Demo only.
#   target/target_chr22.{bed,bim,fam}   REAL public 1000G EUR chr22 subset
#       (100 samples) standing in for the controlled-access validation cohort.
#
# REQUIREMENTS
#   PRS-CS   : git clone https://github.com/getian107/PRScs
#   PRS-CSx  : git clone https://github.com/getian107/PRScsx
#   Python 3 with scipy, numpy, h5py  (h5py must match your numpy — a clean
#       env avoids ABI errors, e.g.  conda create -n prscs -c conda-forge \
#       python=3.10 'numpy<2' scipy h5py)
#   PLINK 1.9 (v1.9.0-b.7.7)
# Point PRSCS_DIR / PRSCSX_DIR at the clones; override PYTHON / PLINK if needed.
# =============================================================================

set -e -o pipefail

PYTHON="${PYTHON:-python}"
PLINK="${PLINK:-plink}"
PRSCS_DIR="${PRSCS_DIR:-$HOME/PRScs}"
PRSCSX_DIR="${PRSCSX_DIR:-$HOME/PRScsx}"

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SD="${SCRIPT_DIR}/sample_data"
OUT="${SCRIPT_DIR}/demo_run"
rm -rf "${OUT}"; mkdir -p "${OUT}"
N_GWAS=50000

echo "=========================================================================="
echo " PRS-CS / PRS-CSx construction — demo smoke test (chr22, 250 SNPs)"
echo "   python     : $(command -v "${PYTHON}" || echo NOT FOUND)"
echo "   plink      : $(command -v "${PLINK}" || echo NOT FOUND)"
echo "   PRSCS_DIR  : ${PRSCS_DIR}"
echo "   PRSCSX_DIR : ${PRSCSX_DIR}"
echo "=========================================================================="
[ -f "${PRSCS_DIR}/PRScs.py" ]   || { echo "ERROR: PRScs.py not under PRSCS_DIR (git clone getian107/PRScs)."; exit 1; }
[ -f "${PRSCSX_DIR}/PRScsx.py" ] || { echo "ERROR: PRScsx.py not under PRSCSX_DIR (git clone getian107/PRScsx)."; exit 1; }

# --- Step 01a: PRS-CS ------------------------------------------------------
echo; echo ">>> Step 01a — PRS-CS (phi=1e-2, chr22)"
"${PYTHON}" "${PRSCS_DIR}/PRScs.py" \
    --ref_dir="${SD}/reference/ldblk_1kg_eur" \
    --bim_prefix="${SD}/target/target_chr22" \
    --sst_file="${SD}/demo_sumstats_chr22.txt" \
    --n_gwas=${N_GWAS} --chrom=22 --phi=1e-2 --seed=9500 \
    --out_dir="${OUT}/prscs" 2>&1 | tail -2
PRSCS_PST="${OUT}/prscs_pst_eff_a1_b0.5_phi1e-02_chr22.txt"

# --- Step 01b: PRS-CSx (EUR) ----------------------------------------------
echo; echo ">>> Step 01b — PRS-CSx (pop=EUR, phi=1e-2, chr22)"
"${PYTHON}" "${PRSCSX_DIR}/PRScsx.py" \
    --ref_dir="${SD}/reference" \
    --bim_prefix="${SD}/target/target_chr22" \
    --sst_file="${SD}/demo_sumstats_chr22.txt" \
    --n_gwas=${N_GWAS} --pop=EUR --chrom=22 --phi=1e-2 --seed=9500 \
    --out_dir="${OUT}" --out_name=prscsx 2>&1 | tail -2
PRSCSX_PST="${OUT}/prscsx_EUR_pst_eff_a1_b0.5_phi1e-02_chr22.txt"

# --- Step 02: merge across chromosomes (single chr here) -------------------
cat "${PRSCS_PST}"  > "${OUT}/prscs_allchr.txt"
cat "${PRSCSX_PST}" > "${OUT}/prscsx_allchr.txt"

# --- Step 03: plink --score ------------------------------------------------
echo; echo ">>> Step 03 — plink --score (SNP=col2, A1=col4, BETA=col6)"
"${PLINK}" --bfile "${SD}/target/target_chr22" \
    --score "${OUT}/prscs_allchr.txt" 2 4 6 sum center \
    --out "${OUT}/prscs_score" 2>&1 | tail -1
"${PLINK}" --bfile "${SD}/target/target_chr22" \
    --score "${OUT}/prscsx_allchr.txt" 2 4 6 sum center \
    --out "${OUT}/prscsx_score" 2>&1 | tail -1

echo; echo ">>> PRS distribution (PRS-CS):"
awk 'NR>1{print $NF}' "${OUT}/prscs_score.profile" | sort -n | \
    awk '{a[NR]=$1} END{print "    n="NR"  min="a[1]"  median="a[int(NR/2)]"  max="a[NR]}'

echo; echo ">>> comparing against the committed reference run ..."
cmpcol(){ paste <(awk '{print $6}' "$1") <(awk '{print $6}' "$2") | \
    awk '{d=($1-$2); if(d<0)d=-d; if(d>m)m=d} END{print m+0}'; }
for pair in "prscs_pst_eff_a1_b0.5_phi1e-02_chr22.txt:${PRSCS_PST}" \
            "prscsx_EUR_pst_eff_a1_b0.5_phi1e-02_chr22.txt:${PRSCSX_PST}"; do
    exp="${SD}/expected_outputs/${pair%%:*}"; got="${pair##*:}"
    if [ -f "$exp" ]; then d=$(cmpcol "$got" "$exp");
        awk -v d="$d" -v f="${pair%%:*}" 'BEGIN{printf "    %-46s max |Δbeta| = %.2e  %s\n", f, d, (d<1e-6?"OK":"DIFF")}'; fi
done

echo; echo "Demo complete. Outputs in ${OUT}/ ."
