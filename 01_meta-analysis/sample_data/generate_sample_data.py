"""
Sample Data Generator
=====================
Produces the small committed `sample_data/` files used by `run_test.py`
when the full GBMI + TAGC source files are not available locally.

Strategy (matches run_test.py with seed=42):
  - 30 common SNPs (in both sample inputs)
  - 70 GBMI-only SNPs
  - 20 TAGC-only SNPs
  -> 100 sample_gbmi rows, 50 sample_tagc rows

After sampling, the IVW meta-analysis + Stage-3 outer join + Stage-4
formatting are run in pure Python (same logic as 02_meta_common.R /
03_merge_all_snps.py / 04_format_downstream.py) and the expected
outputs are written to `sample_data/expected_outputs/`.

Run once whenever you need to regenerate; the produced files are
committed to the repo.

Usage:
  cd 01_meta-analysis/sample_data
  python generate_sample_data.py
"""

import gzip
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

np.random.seed(42)

# --- Source paths (full sumstats live outside the repo). Set the
#     DATA_ROOT env var to the parent directory that contains GBMI/ and TAGC/. ---
DATA_ROOT = os.environ.get("DATA_ROOT")
if not DATA_ROOT:
    raise SystemExit(
        "Set DATA_ROOT to the directory containing GBMI/sumstats/... and "
        "TAGC/TAGC meta-analyses results for asthma risk/... before running.\n"
        "Example:\n  DATA_ROOT=/path/to/hartwell python generate_sample_data.py"
    )
GBMI_PATH = os.path.join(
    DATA_ROOT,
    "GBMI/sumstats/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz",
)
TAGC_PATH = os.path.join(
    DATA_ROOT,
    "TAGC/TAGC meta-analyses results for asthma risk/"
    "TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv",
)

# --- Output locations (inside the repo, committed) ---
SAMPLE_DIR = os.path.dirname(os.path.abspath(__file__))
EXPECTED_DIR = os.path.join(SAMPLE_DIR, "expected_outputs")
os.makedirs(EXPECTED_DIR, exist_ok=True)

SAMPLE_GBMI = os.path.join(SAMPLE_DIR, "sample_gbmi.txt.gz")
SAMPLE_TAGC = os.path.join(SAMPLE_DIR, "sample_tagc.tsv")
EXP_COMMON = os.path.join(EXPECTED_DIR, "meta_analysis_results_common.csv")
EXP_OUTPUT = os.path.join(EXPECTED_DIR, "meta_analysis_output.csv")
EXP_COMBINE_STAT = os.path.join(EXPECTED_DIR, "meta_analysis_combineStat.txt")
EXP_FORMATTED = os.path.join(EXPECTED_DIR, "formatted_meta_analysis.txt")
EXP_P_COMPLETE = os.path.join(EXPECTED_DIR, "meta_analysis_p_complete.txt.gz")


# ============================================================
# STEP 1: Sample SNPs from the full sumstats
# ============================================================
print("Loading GBMI EUR (this may take a moment)...")
with gzip.open(GBMI_PATH, "rt") as f:
    gbmi_full = pd.read_csv(f, sep="\t")
print(f"  GBMI loaded: {len(gbmi_full):,} rows")

print("Loading TAGC...")
tagc_full = pd.read_csv(TAGC_PATH, sep="\t")
print(f"  TAGC loaded: {len(tagc_full):,} rows")

gbmi_rsids = set(gbmi_full["rsid"].dropna())
tagc_rsids = set(tagc_full["rsid"].dropna())
common_rsids = gbmi_rsids & tagc_rsids
print(f"  Common rsids: {len(common_rsids):,}")

common_sample = np.random.choice(list(common_rsids), size=30, replace=False)
gbmi_only_rsids = gbmi_rsids - tagc_rsids
gbmi_only_sample = np.random.choice(list(gbmi_only_rsids), size=70, replace=False)
tagc_only_rsids = tagc_rsids - gbmi_rsids
tagc_only_sample = np.random.choice(list(tagc_only_rsids), size=20, replace=False)

gbmi_snps = list(common_sample) + list(gbmi_only_sample)
sample_gbmi = (
    gbmi_full[gbmi_full["rsid"].isin(gbmi_snps)]
    .drop_duplicates(subset="rsid")
    .head(100)
    .copy()
)
tagc_snps = list(common_sample) + list(tagc_only_sample)
sample_tagc = (
    tagc_full[tagc_full["rsid"].isin(tagc_snps)]
    .drop_duplicates(subset="rsid")
    .head(50)
    .copy()
)

sample_gbmi.to_csv(SAMPLE_GBMI, sep="\t", index=False, compression="gzip")
sample_tagc.to_csv(SAMPLE_TAGC, sep="\t", index=False)

del gbmi_full, tagc_full

actual_common = set(sample_gbmi["rsid"]) & set(sample_tagc["rsid"])
print(
    f"\nSample inputs written:"
    f"\n  sample_gbmi.txt.gz: {len(sample_gbmi)} rows"
    f"\n  sample_tagc.tsv:    {len(sample_tagc)} rows"
    f"\n  Overlapping rsids:  {len(actual_common)}"
)


# ============================================================
# STEP 2: Stage 2 -- IVW fixed-effects meta-analysis (common SNPs)
# ============================================================
snps = sample_tagc[
    ["chr", "rsid", "reference_allele", "alternate_allele",
     "European_ancestry_beta_fix", "European_ancestry_se_fix"]
].merge(
    sample_gbmi[["#CHR", "rsid", "REF", "ALT", "inv_var_meta_beta", "inv_var_meta_sebeta"]],
    left_on=["chr", "rsid"], right_on=["#CHR", "rsid"], how="inner",
)

concordant = (snps["alternate_allele"] == snps["ALT"]) & (snps["reference_allele"] == snps["REF"])
flipped = (snps["alternate_allele"] == snps["REF"]) & (snps["reference_allele"] == snps["ALT"])
mismatch = ~concordant & ~flipped
snps.loc[flipped, "inv_var_meta_beta"] = -snps.loc[flipped, "inv_var_meta_beta"]
snps = snps[~mismatch].copy()
snps = snps.drop(columns=["reference_allele", "alternate_allele", "REF", "ALT"]).rename(columns={
    "European_ancestry_beta_fix": "beta_tagc",
    "European_ancestry_se_fix": "se_tagc",
    "inv_var_meta_beta": "beta_gbmi",
    "inv_var_meta_sebeta": "se_gbmi",
})

w_tagc = 1.0 / snps["se_tagc"] ** 2
w_gbmi = 1.0 / snps["se_gbmi"] ** 2
snps["combined_beta"] = (w_tagc * snps["beta_tagc"] + w_gbmi * snps["beta_gbmi"]) / (w_tagc + w_gbmi)
snps["combined_se"] = np.sqrt(1.0 / (w_tagc + w_gbmi))

stage2 = snps[["chr", "rsid", "beta_tagc", "se_tagc", "beta_gbmi", "se_gbmi",
               "combined_beta", "combined_se"]]
stage2.to_csv(EXP_COMMON, index=False)
print(f"  meta_analysis_results_common.csv: {len(stage2)} SNPs")


# ============================================================
# STEP 3: Stage 3 -- recover study-specific SNPs (outer join)
# ============================================================
merged = stage2[["chr", "rsid", "combined_beta", "combined_se"]].merge(
    sample_gbmi[["#CHR", "rsid", "inv_var_meta_beta", "inv_var_meta_sebeta"]],
    left_on=["chr", "rsid"], right_on=["#CHR", "rsid"], how="outer",
)
merged["combined_beta"] = merged["combined_beta"].fillna(merged["inv_var_meta_beta"])
merged["combined_se"] = merged["combined_se"].fillna(merged["inv_var_meta_sebeta"])
merged["chr"] = merged["chr"].fillna(merged["#CHR"])

final = merged.merge(
    sample_tagc[["chr", "rsid", "European_ancestry_beta_fix", "European_ancestry_se_fix"]],
    on=["chr", "rsid"], how="outer",
)
final["combined_beta"] = final["combined_beta"].fillna(final["European_ancestry_beta_fix"])
final["combined_se"] = final["combined_se"].fillna(final["European_ancestry_se_fix"])

stage3 = final.drop(columns=["#CHR"]).rename(columns={
    "inv_var_meta_beta": "gbmi_beta",
    "inv_var_meta_sebeta": "gbmi_se",
    "European_ancestry_beta_fix": "tagc_beta",
    "European_ancestry_se_fix": "tagc_se",
})
stage3.to_csv(EXP_OUTPUT, index=False)
print(f"  meta_analysis_output.csv:         {len(stage3)} SNPs")


# ============================================================
# STEP 4: Stage 4 -- format for downstream
# ============================================================
summary_stats = stage3.copy()
summary_stats["Z"] = summary_stats["combined_beta"] / summary_stats["combined_se"]

# Output 1: combineStat.txt (BETA, SE)
combine_stat = (
    summary_stats[["rsid", "combined_beta", "combined_se"]]
    .rename(columns={"rsid": "SNP", "combined_beta": "BETA", "combined_se": "SE"})
)
combine_stat.to_csv(EXP_COMBINE_STAT, sep="\t", index=False)
print(f"  meta_analysis_combineStat.txt:    {len(combine_stat)} SNPs")

# Output 2: formatted_meta_analysis.txt (SNP, Z, A2, A1)
formatted = (
    summary_stats[["rsid", "combined_beta", "combined_se", "Z"]]
    .rename(columns={"rsid": "SNP", "combined_beta": "BETA", "combined_se": "SE"})
)
formatted = formatted.merge(
    sample_tagc[["rsid", "reference_allele", "alternate_allele"]],
    left_on="SNP", right_on="rsid", how="left",
)
formatted["A1"] = formatted["alternate_allele"]
formatted["A2"] = formatted["reference_allele"]
formatted.drop(columns=["rsid", "reference_allele", "alternate_allele"], inplace=True)

missing_mask = formatted["A1"].isna()
ref_dict = sample_gbmi.set_index("rsid")["REF"].to_dict()
alt_dict = sample_gbmi.set_index("rsid")["ALT"].to_dict()
formatted.loc[missing_mask, "A1"] = formatted.loc[missing_mask, "SNP"].map(alt_dict)
formatted.loc[missing_mask, "A2"] = formatted.loc[missing_mask, "SNP"].map(ref_dict)
formatted = formatted.dropna(subset=["A1", "A2"])

formatted[["SNP", "Z", "A2", "A1"]].to_csv(EXP_FORMATTED, sep="\t", index=False)
print(f"  formatted_meta_analysis.txt:      {len(formatted)} SNPs")

# Output 3: meta_analysis_p_complete.txt.gz (with positions and p-values)
p_complete = formatted.copy()
p_complete = p_complete.merge(
    sample_gbmi[["rsid", "#CHR", "POS"]], left_on="SNP", right_on="rsid", how="left"
).drop(columns=["rsid"])
p_complete = p_complete.merge(
    sample_tagc[["rsid", "chr", "position"]], left_on="SNP", right_on="rsid", how="left"
)
p_complete["#CHR"] = p_complete["#CHR"].fillna(p_complete["chr"])
p_complete["POS"] = p_complete["POS"].fillna(p_complete["position"])
p_complete.drop(columns=["rsid", "chr", "position"], inplace=True)

p_complete["p_value"] = 2 * norm.sf(np.abs(p_complete["Z"]))
p_complete["-log10(P-Value)"] = -np.log10(p_complete["p_value"])
p_complete = p_complete.dropna(subset=["#CHR", "POS"])
p_complete = p_complete[
    ["SNP", "A1", "A2", "Z", "#CHR", "POS", "BETA", "SE", "p_value", "-log10(P-Value)"]
]
p_complete.to_csv(EXP_P_COMPLETE, sep="\t", index=False, compression="gzip")
print(f"  meta_analysis_p_complete.txt.gz:  {len(p_complete)} SNPs")

print("\nDone. Commit the contents of sample_data/ to the repo.")
