"""
Stage 3: Recover study-specific SNPs via outer joins

Input:
  - output/meta_analysis_results_common.csv  (from Stage 2)
  - GBMI EUR raw sumstats
  - TAGC raw sumstats

Output:
  - output/meta_analysis_output.csv
    Columns: chr, rsid, combined_beta, combined_se,
             gbmi_beta, gbmi_se, tagc_beta, tagc_se

For SNPs in both studies:  combined = IVW meta-analysis estimate
For GBMI-only SNPs:        combined = GBMI estimate
For TAGC-only SNPs:        combined = TAGC estimate
"""

import pandas as pd
import gzip

# --- Configuration ---
COMMON_PATH = "output/meta_analysis_results_common.csv"
GBMI_PATH = "../GBMI/sumstats/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
TAGC_PATH = "../TAGC/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv"
OUTPUT_PATH = "output/meta_analysis_output.csv"

# --- Load data ---
print("Loading data...")
common = pd.read_csv(COMMON_PATH)

with gzip.open(GBMI_PATH, "rt") as f:
    eur = pd.read_csv(f, sep="\t")

tagc = pd.read_csv(TAGC_PATH, sep="\t")

print(f"  Common SNPs (TAGC+GBMI): {len(common)}")
print(f"  GBMI EUR SNPs:           {len(eur)}")
print(f"  TAGC SNPs:               {len(tagc)}")

# --- Step 1: Outer join with GBMI EUR to recover GBMI-only SNPs ---
merged = common[["chr", "rsid", "combined_beta", "combined_se"]].merge(
    eur[["#CHR", "rsid", "inv_var_meta_beta", "inv_var_meta_sebeta"]],
    left_on=["chr", "rsid"],
    right_on=["#CHR", "rsid"],
    how="outer",
)

# Fill missing combined values with GBMI-only estimates
merged["combined_beta"] = merged["combined_beta"].fillna(merged["inv_var_meta_beta"])
merged["combined_se"] = merged["combined_se"].fillna(merged["inv_var_meta_sebeta"])

# Fill chr from #CHR for GBMI-only SNPs
merged["chr"] = merged["chr"].fillna(merged["#CHR"])

# --- Step 2: Outer join with TAGC to recover TAGC-only SNPs ---
final = merged.merge(
    tagc[["chr", "rsid", "European_ancestry_beta_fix", "European_ancestry_se_fix"]],
    on=["chr", "rsid"],
    how="outer",
)

# Fill missing combined values with TAGC-only estimates
final["combined_beta"] = final["combined_beta"].fillna(
    final["European_ancestry_beta_fix"]
)
final["combined_se"] = final["combined_se"].fillna(
    final["European_ancestry_se_fix"]
)

# --- Clean up and save ---
output = final.drop(columns=["#CHR"]).rename(
    columns={
        "inv_var_meta_beta": "gbmi_beta",
        "inv_var_meta_sebeta": "gbmi_se",
        "European_ancestry_beta_fix": "tagc_beta",
        "European_ancestry_se_fix": "tagc_se",
    }
)

output.to_csv(OUTPUT_PATH, index=False)
print(f"\nSaved {len(output)} total SNPs to {OUTPUT_PATH}")
print(f"  Both studies:  {common.shape[0]}")
print(f"  GBMI-only:     {output['tagc_beta'].isna().sum() - (len(output) - len(eur) - len(tagc) + common.shape[0])}")
print(f"  TAGC-only:     {output['gbmi_beta'].isna().sum() - (len(output) - len(eur) - len(tagc) + common.shape[0])}")
