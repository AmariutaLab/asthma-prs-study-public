"""
Pipeline Test Script
====================
Runs the meta-analysis pipeline (Stages 2-4) and verifies the outputs.

Two modes, auto-selected based on what's available locally:

  REAL mode (full sumstats present)
    1. Sample 100 GBMI + 50 TAGC SNPs from the full sumstats (seed=42)
    2. Run the pipeline on the sample
    3. Compare against the real full pipeline outputs

  SAMPLE mode (full sumstats absent)
    1. Load the committed sample inputs from sample_data/
    2. Run the pipeline on them
    3. Compare against the committed expected outputs in sample_data/expected_outputs/

SAMPLE mode lets anyone clone the repo and run the test without needing
the ~1 GB GBMI/TAGC source files. The two modes share Stages 2-4 verbatim.

Usage: python run_test.py
"""

import pandas as pd
import gzip
import numpy as np
from scipy.stats import norm
import os

np.random.seed(42)

# --- Paths ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
HARTWELL = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
GBMI_PATH = os.path.join(HARTWELL, "GBMI/sumstats/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz")
TAGC_PATH = os.path.join(HARTWELL, "TAGC/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv")

# Real results for comparison (REAL mode only)
REAL_COMMON = os.path.join(HARTWELL, "meta_analysis/meta_analysis_results_common.csv")
REAL_OUTPUT = os.path.join(HARTWELL, "TWAS/meta_analysis_output.csv")
REAL_COMBINE_STAT = os.path.join(HARTWELL, "TWAS/meta_analysis_combineStat.txt")
REAL_FORMATTED = os.path.join(HARTWELL, "TWAS/formatted_meta_analysis.txt")
REAL_P_COMPLETE = os.path.join(HARTWELL, "Clumping/meta_analysis_p_complete.txt.gz")

# Committed sample inputs + expected outputs (SAMPLE mode only)
SAMPLE_DIR = os.path.join(SCRIPT_DIR, "sample_data")
SAMPLE_GBMI = os.path.join(SAMPLE_DIR, "sample_gbmi.txt.gz")
SAMPLE_TAGC = os.path.join(SAMPLE_DIR, "sample_tagc.tsv")
SAMPLE_EXPECTED_DIR = os.path.join(SAMPLE_DIR, "expected_outputs")
SAMPLE_COMMON = os.path.join(SAMPLE_EXPECTED_DIR, "meta_analysis_results_common.csv")
SAMPLE_OUTPUT = os.path.join(SAMPLE_EXPECTED_DIR, "meta_analysis_output.csv")
SAMPLE_COMBINE_STAT = os.path.join(SAMPLE_EXPECTED_DIR, "meta_analysis_combineStat.txt")
SAMPLE_FORMATTED = os.path.join(SAMPLE_EXPECTED_DIR, "formatted_meta_analysis.txt")
SAMPLE_P_COMPLETE = os.path.join(SAMPLE_EXPECTED_DIR, "meta_analysis_p_complete.txt.gz")

# --- Mode selection ---
USE_REAL = os.path.exists(GBMI_PATH) and os.path.exists(TAGC_PATH)
MODE = "REAL" if USE_REAL else "SAMPLE"

# Per-mode comparison targets
if USE_REAL:
    EXPECTED_COMMON = REAL_COMMON
    EXPECTED_OUTPUT = REAL_OUTPUT
    EXPECTED_COMBINE_STAT = REAL_COMBINE_STAT
    EXPECTED_FORMATTED = REAL_FORMATTED
    EXPECTED_P_COMPLETE = REAL_P_COMPLETE
else:
    EXPECTED_COMMON = SAMPLE_COMMON
    EXPECTED_OUTPUT = SAMPLE_OUTPUT
    EXPECTED_COMBINE_STAT = SAMPLE_COMBINE_STAT
    EXPECTED_FORMATTED = SAMPLE_FORMATTED
    EXPECTED_P_COMPLETE = SAMPLE_P_COMPLETE

# Test output directory
TEST_DIR = os.path.join(SCRIPT_DIR, "test")
os.makedirs(TEST_DIR, exist_ok=True)

# ============================================================
# STEP 1: Acquire test inputs
# ============================================================
print("=" * 60)
print(f"STEP 1: Acquire test inputs  [mode = {MODE}]")
print("=" * 60)

test_gbmi_path = os.path.join(TEST_DIR, "test_gbmi.txt.gz")
test_tagc_path = os.path.join(TEST_DIR, "test_tagc.tsv")

if USE_REAL:
    print("Loading GBMI EUR (this may take a moment)...")
    with gzip.open(GBMI_PATH, "rt") as f:
        gbmi_full = pd.read_csv(f, sep="\t")
    print(f"  GBMI loaded: {len(gbmi_full)} rows")

    print("Loading TAGC...")
    tagc_full = pd.read_csv(TAGC_PATH, sep="\t")
    print(f"  TAGC loaded: {len(tagc_full)} rows")

    # Find common rsids (non-null)
    gbmi_rsids = set(gbmi_full["rsid"].dropna())
    tagc_rsids = set(tagc_full["rsid"].dropna())
    common_rsids = gbmi_rsids & tagc_rsids
    print(f"  Common rsids: {len(common_rsids)}")

    # Sample with guaranteed overlap:
    #   30 common SNPs (in both test files)
    #   70 GBMI-only SNPs
    #   20 TAGC-only SNPs
    # Total: 100 GBMI rows, 50 TAGC rows
    common_sample = np.random.choice(list(common_rsids), size=30, replace=False)
    gbmi_only_rsids = gbmi_rsids - tagc_rsids
    gbmi_only_sample = np.random.choice(list(gbmi_only_rsids), size=70, replace=False)
    tagc_only_rsids = tagc_rsids - gbmi_rsids
    tagc_only_sample = np.random.choice(list(tagc_only_rsids), size=20, replace=False)

    gbmi_snps = list(common_sample) + list(gbmi_only_sample)
    test_gbmi = gbmi_full[gbmi_full["rsid"].isin(gbmi_snps)].copy()
    test_gbmi = test_gbmi.drop_duplicates(subset="rsid").head(100)

    tagc_snps = list(common_sample) + list(tagc_only_sample)
    test_tagc = tagc_full[tagc_full["rsid"].isin(tagc_snps)].copy()
    test_tagc = test_tagc.drop_duplicates(subset="rsid").head(50)

    test_gbmi.to_csv(test_gbmi_path, sep="\t", index=False, compression="gzip")
    test_tagc.to_csv(test_tagc_path, sep="\t", index=False)

    del gbmi_full
else:
    print("Full GBMI / TAGC sumstats not found; using committed sample_data/.")
    print(f"  GBMI source: {SAMPLE_GBMI}")
    print(f"  TAGC source: {SAMPLE_TAGC}")
    with gzip.open(SAMPLE_GBMI, "rt") as f:
        test_gbmi = pd.read_csv(f, sep="\t")
    test_tagc = pd.read_csv(SAMPLE_TAGC, sep="\t")
    test_gbmi.to_csv(test_gbmi_path, sep="\t", index=False, compression="gzip")
    test_tagc.to_csv(test_tagc_path, sep="\t", index=False)

actual_common = set(test_gbmi["rsid"]) & set(test_tagc["rsid"])
print(f"\nTest data ready:")
print(f"  test_gbmi.txt.gz:  {len(test_gbmi)} rows")
print(f"  test_tagc.tsv:     {len(test_tagc)} rows")
print(f"  Overlapping rsids: {len(actual_common)}")

# ============================================================
# STEP 2: Run Stage 2 - Meta-analysis of common SNPs
# (Reimplemented in Python; identical IVW fixed-effects math)
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Stage 2 - Meta-analysis of common SNPs")
print("=" * 60)

# Reload test data
with gzip.open(test_gbmi_path, "rt") as f:
    t_gbmi = pd.read_csv(f, sep="\t")
t_tagc = pd.read_csv(test_tagc_path, sep="\t")

# Inner join on common SNPs (same as 02_meta_common.R)
snps = t_tagc[["chr", "rsid", "reference_allele", "alternate_allele",
               "European_ancestry_beta_fix", "European_ancestry_se_fix"]].merge(
    t_gbmi[["#CHR", "rsid", "REF", "ALT", "inv_var_meta_beta", "inv_var_meta_sebeta"]],
    left_on=["chr", "rsid"],
    right_on=["#CHR", "rsid"],
    how="inner",
)

# Allele concordance check (same as 02_meta_common.R)
concordant = (snps["alternate_allele"] == snps["ALT"]) & (snps["reference_allele"] == snps["REF"])
flipped = (snps["alternate_allele"] == snps["REF"]) & (snps["reference_allele"] == snps["ALT"])
mismatch = ~concordant & ~flipped

print(f"  Allele check: {concordant.sum()} concordant, {flipped.sum()} flipped, {mismatch.sum()} mismatch")

# Flip GBMI beta for reversed alleles, remove mismatches
snps.loc[flipped, "inv_var_meta_beta"] = -snps.loc[flipped, "inv_var_meta_beta"]
snps = snps[~mismatch].copy()

snps = snps.drop(columns=["reference_allele", "alternate_allele", "REF", "ALT"]).rename(columns={
    "European_ancestry_beta_fix": "beta_tagc",
    "European_ancestry_se_fix": "se_tagc",
    "inv_var_meta_beta": "beta_gbmi",
    "inv_var_meta_sebeta": "se_gbmi",
})

# IVW fixed-effects meta-analysis (same formula as metafor::rma with method="FE")
w_tagc = 1.0 / snps["se_tagc"] ** 2
w_gbmi = 1.0 / snps["se_gbmi"] ** 2
snps["combined_beta"] = (w_tagc * snps["beta_tagc"] + w_gbmi * snps["beta_gbmi"]) / (w_tagc + w_gbmi)
snps["combined_se"] = np.sqrt(1.0 / (w_tagc + w_gbmi))

stage2_out = snps[["chr", "rsid", "beta_tagc", "se_tagc", "beta_gbmi", "se_gbmi", "combined_beta", "combined_se"]]
stage2_path = os.path.join(TEST_DIR, "meta_analysis_results_common.csv")
stage2_out.to_csv(stage2_path, index=False)
print(f"  Common SNPs meta-analyzed: {len(stage2_out)}")

# ============================================================
# STEP 3: Stage 3 - Merge study-specific SNPs
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Stage 3 - Recover study-specific SNPs")
print("=" * 60)

common_df = pd.read_csv(stage2_path)

# Outer join with GBMI
merged = common_df[["chr", "rsid", "combined_beta", "combined_se"]].merge(
    t_gbmi[["#CHR", "rsid", "inv_var_meta_beta", "inv_var_meta_sebeta"]],
    left_on=["chr", "rsid"],
    right_on=["#CHR", "rsid"],
    how="outer",
)
merged["combined_beta"] = merged["combined_beta"].fillna(merged["inv_var_meta_beta"])
merged["combined_se"] = merged["combined_se"].fillna(merged["inv_var_meta_sebeta"])
merged["chr"] = merged["chr"].fillna(merged["#CHR"])

# Outer join with TAGC
final = merged.merge(
    t_tagc[["chr", "rsid", "European_ancestry_beta_fix", "European_ancestry_se_fix"]],
    on=["chr", "rsid"],
    how="outer",
)
final["combined_beta"] = final["combined_beta"].fillna(final["European_ancestry_beta_fix"])
final["combined_se"] = final["combined_se"].fillna(final["European_ancestry_se_fix"])

output = final.drop(columns=["#CHR"]).rename(columns={
    "inv_var_meta_beta": "gbmi_beta",
    "inv_var_meta_sebeta": "gbmi_se",
    "European_ancestry_beta_fix": "tagc_beta",
    "European_ancestry_se_fix": "tagc_se",
})

stage3_path = os.path.join(TEST_DIR, "meta_analysis_output.csv")
output.to_csv(stage3_path, index=False)

n_both = common_df.shape[0]
n_gbmi_only = output["tagc_beta"].isna().sum()
n_tagc_only = output["gbmi_beta"].isna().sum()
print(f"  Total SNPs:    {len(output)}")
print(f"  Both studies:  {n_both}")
print(f"  GBMI-only:     {n_gbmi_only}")
print(f"  TAGC-only:     {n_tagc_only}")

# ============================================================
# STEP 4: Stage 4 - Format for downstream
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Stage 4 - Format for downstream analyses")
print("=" * 60)

summary_stats = pd.read_csv(stage3_path, low_memory=False)
summary_stats["Z"] = summary_stats["combined_beta"] / summary_stats["combined_se"]

# Output 1: meta_analysis_combineStat.txt
combine_stat = summary_stats[["rsid", "combined_beta", "combined_se"]].copy()
combine_stat.rename(columns={"rsid": "SNP", "combined_beta": "BETA", "combined_se": "SE"}, inplace=True)
combine_stat_path = os.path.join(TEST_DIR, "meta_analysis_combineStat.txt")
combine_stat.to_csv(combine_stat_path, sep="\t", index=False)
print(f"  meta_analysis_combineStat.txt: {len(combine_stat)} SNPs")

# Assign alleles
formatted = summary_stats[["rsid", "combined_beta", "combined_se", "Z"]].copy()
formatted.rename(columns={"rsid": "SNP", "combined_beta": "BETA", "combined_se": "SE"}, inplace=True)

# From TAGC
formatted = formatted.merge(
    t_tagc[["rsid", "reference_allele", "alternate_allele"]],
    left_on="SNP", right_on="rsid", how="left",
)
formatted["A1"] = formatted["alternate_allele"]  # ALT = effect allele
formatted["A2"] = formatted["reference_allele"]  # REF
formatted.drop(columns=["rsid", "reference_allele", "alternate_allele"], inplace=True)

# From GBMI (fill remaining)
missing_mask = formatted["A1"].isna()
ref_dict = t_gbmi.set_index("rsid")["REF"].to_dict()
alt_dict = t_gbmi.set_index("rsid")["ALT"].to_dict()
formatted.loc[missing_mask, "A1"] = formatted.loc[missing_mask, "SNP"].map(alt_dict)
formatted.loc[missing_mask, "A2"] = formatted.loc[missing_mask, "SNP"].map(ref_dict)
formatted = formatted.dropna(subset=["A1", "A2"])

# Output 2: formatted_meta_analysis.txt
formatted_path = os.path.join(TEST_DIR, "formatted_meta_analysis.txt")
formatted[["SNP", "Z", "A2", "A1"]].to_csv(formatted_path, sep="\t", index=False)
print(f"  formatted_meta_analysis.txt:   {len(formatted)} SNPs")

# Output 3: meta_analysis_p_complete.txt.gz
p_complete = formatted.copy()

# Positions from GBMI
p_complete = p_complete.merge(
    t_gbmi[["rsid", "#CHR", "POS"]], left_on="SNP", right_on="rsid", how="left"
)
p_complete.drop(columns=["rsid"], inplace=True)

# Fill from TAGC
p_complete = p_complete.merge(
    t_tagc[["rsid", "chr", "position"]], left_on="SNP", right_on="rsid", how="left"
)
p_complete["#CHR"] = p_complete["#CHR"].fillna(p_complete["chr"])
p_complete["POS"] = p_complete["POS"].fillna(p_complete["position"])
p_complete.drop(columns=["rsid", "chr", "position"], inplace=True)

p_complete["p_value"] = 2 * norm.sf(np.abs(p_complete["Z"]))
p_complete["-log10(P-Value)"] = -np.log10(p_complete["p_value"])
p_complete = p_complete.dropna(subset=["#CHR", "POS"])
p_complete = p_complete[["SNP", "A1", "A2", "Z", "#CHR", "POS", "BETA", "SE", "p_value", "-log10(P-Value)"]]

p_complete_path = os.path.join(TEST_DIR, "meta_analysis_p_complete.txt.gz")
p_complete.to_csv(p_complete_path, sep="\t", index=False, compression="gzip")
print(f"  meta_analysis_p_complete.txt.gz: {len(p_complete)} SNPs")

# ============================================================
# STEP 5: Compare with real results
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: Comparing test outputs with real results")
print("=" * 60)

all_pass = True

def compare_df(test_df, real_df, key_col, value_cols, label):
    """Compare matching rows between test and real DataFrames."""
    global all_pass
    merged = test_df.merge(real_df, on=key_col, suffixes=("_test", "_real"))
    if len(merged) == 0:
        print(f"\n  [{label}] WARNING: No matching SNPs found!")
        all_pass = False
        return

    print(f"\n  [{label}] Matched {len(merged)} / {len(test_df)} SNPs")

    for col in value_cols:
        tc = f"{col}_test"
        rc = f"{col}_real"
        if tc not in merged.columns or rc not in merged.columns:
            # Columns might not have suffix if key_col == col
            continue
        # Skip non-numeric
        if not pd.api.types.is_numeric_dtype(merged[tc]):
            match = (merged[tc] == merged[rc]).all()
            print(f"    {col}: {'MATCH' if match else 'MISMATCH'}")
            if not match:
                all_pass = False
                mismatches = merged[merged[tc] != merged[rc]][[key_col, tc, rc]].head(3)
                print(f"      First mismatches:\n{mismatches.to_string(index=False)}")
            continue

        # The REAL_P_COMPLETE file was generated with a one-sided p_value
        # formula; the pipeline (and the SAMPLE expected outputs) use two-sided,
        # so in REAL mode test_p_value should equal 2 * real_p_value.
        if col == "p_value" and USE_REAL:
            diff = (merged[tc] - 2 * merged[rc]).abs()
        else:
            diff = (merged[tc] - merged[rc]).abs()
        max_diff = diff.max()
        tol = 1e-6
        status = "PASS" if max_diff < tol else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"    {col}: max_diff={max_diff:.2e}  {status}")
        if status == "FAIL":
            worst = merged.loc[diff.idxmax(), [key_col, tc, rc]]
            print(f"      Worst: {worst.to_dict()}")


# --- Compare Stage 2: meta_analysis_results_common.csv ---
print("\n--- Stage 2: meta_analysis_results_common.csv ---")
test_common = pd.read_csv(stage2_path)
real_common = pd.read_csv(EXPECTED_COMMON)
compare_df(test_common, real_common, "rsid",
           ["combined_beta", "combined_se", "beta_tagc", "se_tagc", "beta_gbmi", "se_gbmi"],
           "Stage 2")


# --- Compare Stage 3: meta_analysis_output.csv ---
print("\n--- Stage 3: meta_analysis_output.csv ---")
test_output = pd.read_csv(stage3_path, low_memory=False)
real_output = pd.read_csv(EXPECTED_OUTPUT, low_memory=False)
compare_df(test_output, real_output, "rsid",
           ["combined_beta", "combined_se", "gbmi_beta", "gbmi_se", "tagc_beta", "tagc_se"],
           "Stage 3")


# --- Compare Stage 4a: meta_analysis_combineStat.txt ---
print("\n--- Stage 4a: meta_analysis_combineStat.txt ---")
test_cs = pd.read_csv(combine_stat_path, sep="\t")
real_cs = pd.read_csv(EXPECTED_COMBINE_STAT, sep="\t", low_memory=False)
compare_df(test_cs, real_cs, "SNP", ["BETA", "SE"], "Stage 4a")


# --- Compare Stage 4b: formatted_meta_analysis.txt ---
print("\n--- Stage 4b: formatted_meta_analysis.txt ---")
test_fmt = pd.read_csv(formatted_path, sep="\t")
real_fmt = pd.read_csv(EXPECTED_FORMATTED, sep="\t")
compare_df(test_fmt, real_fmt, "SNP", ["Z", "A1", "A2"], "Stage 4b")


# --- Compare Stage 4c: meta_analysis_p_complete.txt.gz ---
print("\n--- Stage 4c: meta_analysis_p_complete.txt.gz ---")
test_pc = pd.read_csv(p_complete_path, sep="\t")
with gzip.open(EXPECTED_P_COMPLETE, "rt") as f:
    real_pc = pd.read_csv(f, sep="\t")
compare_df(test_pc, real_pc, "SNP",
           ["A1", "A2", "Z", "#CHR", "POS", "BETA", "SE", "p_value"],
           "Stage 4c")


# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
if all_pass:
    print("ALL COMPARISONS PASSED")
else:
    print("SOME COMPARISONS FAILED - check details above")
print("=" * 60)
