"""
Stage 4: Format meta-analysis results for downstream analyses (TWAS, PRS)

Input:
  - output/meta_analysis_output.csv  (from Stage 3)
  - GBMI EUR raw sumstats            (for alleles and genomic positions)
  - TAGC raw sumstats                (for alleles and genomic positions)

Output:
  - output/meta_analysis_combineStat.txt   (SNP, BETA, SE)
  - output/formatted_meta_analysis.txt     (SNP, Z, A2, A1)  for TWAS
  - output/meta_analysis_p_complete.txt.gz (SNP, A1, A2, Z, #CHR, POS,
                                            BETA, SE, p_value, -log10(P-Value))
                                           for Clumping/PRS

Allele convention throughout:
  A1 = ALT = effect allele  (beta is relative to this allele)
  A2 = REF = reference allele
"""

import pandas as pd
import gzip
import numpy as np
from scipy.stats import norm

# --- Configuration ---
META_OUTPUT_PATH = "output/meta_analysis_output.csv"
GBMI_PATH = "../GBMI/sumstats/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
TAGC_PATH = (
    "../TAGC/TAGC meta-analyses results for asthma risk/"
    "TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv"
)

OUTPUT_COMBINE_STAT = "output/meta_analysis_combineStat.txt"
OUTPUT_FORMATTED = "output/formatted_meta_analysis.txt"
OUTPUT_P_COMPLETE = "output/meta_analysis_p_complete.txt.gz"

# ============================================================
# Load data
# ============================================================
print("Loading meta-analysis output...")
summary_stats = pd.read_csv(META_OUTPUT_PATH, low_memory=False)
summary_stats["Z"] = summary_stats["combined_beta"] / summary_stats["combined_se"]

print("Loading TAGC...")
tagc = pd.read_csv(TAGC_PATH, sep="\t")

print("Loading GBMI EUR...")
with gzip.open(GBMI_PATH, "rt") as f:
    eur = pd.read_csv(f, sep="\t")

# ============================================================
# Output 1: meta_analysis_combineStat.txt  (SNP, BETA, SE)
# ============================================================
combine_stat = summary_stats[["rsid", "combined_beta", "combined_se"]].copy()
combine_stat.rename(
    columns={"rsid": "SNP", "combined_beta": "BETA", "combined_se": "SE"},
    inplace=True,
)
combine_stat.to_csv(OUTPUT_COMBINE_STAT, sep="\t", index=False)
print(f"\nSaved {OUTPUT_COMBINE_STAT}: {len(combine_stat)} SNPs")

# ============================================================
# Assign alleles: A1 = ALT (effect), A2 = REF
# ============================================================
formatted = summary_stats[["rsid", "combined_beta", "combined_se", "Z"]].copy()
formatted.rename(
    columns={"rsid": "SNP", "combined_beta": "BETA", "combined_se": "SE"},
    inplace=True,
)

# Step 1: Get alleles from TAGC (alternate_allele -> A1, reference_allele -> A2)
formatted = formatted.merge(
    tagc[["rsid", "reference_allele", "alternate_allele"]],
    left_on="SNP",
    right_on="rsid",
    how="left",
)
formatted["A1"] = formatted["alternate_allele"]  # A1 = ALT = effect allele
formatted["A2"] = formatted["reference_allele"]  # A2 = REF
formatted.drop(
    columns=["rsid", "reference_allele", "alternate_allele"], inplace=True
)

tagc_filled = (~formatted["A1"].isna()).sum()
print(f"\nAllele assignment:")
print(f"  From TAGC:        {tagc_filled}")

# Step 2: Fill remaining alleles from GBMI (ALT -> A1, REF -> A2)
missing_mask = formatted["A1"].isna()
print(f"  Needing GBMI:     {missing_mask.sum()}")

ref_dict = eur.set_index("rsid")["REF"].to_dict()
alt_dict = eur.set_index("rsid")["ALT"].to_dict()

formatted.loc[missing_mask, "A1"] = formatted.loc[missing_mask, "SNP"].map(alt_dict)
formatted.loc[missing_mask, "A2"] = formatted.loc[missing_mask, "SNP"].map(ref_dict)

# Step 3: Drop SNPs still missing alleles
before = len(formatted)
formatted = formatted.dropna(subset=["A1", "A2"])
print(f"  Dropped (no alleles): {before - len(formatted)}")
print(f"  Final SNP count:      {len(formatted)}")

# ============================================================
# Output 2: formatted_meta_analysis.txt  (SNP, Z, A2, A1)
# ============================================================
formatted_out = formatted[["SNP", "Z", "A2", "A1"]]
formatted_out.to_csv(OUTPUT_FORMATTED, sep="\t", index=False)
print(f"\nSaved {OUTPUT_FORMATTED}: {len(formatted_out)} SNPs")

# ============================================================
# Output 3: meta_analysis_p_complete.txt.gz
#   (SNP, A1, A2, Z, #CHR, POS, BETA, SE, p_value, -log10(P-Value))
# ============================================================
print("\nBuilding meta_analysis_p_complete...")

# Start from formatted (has SNP, BETA, SE, Z, A1, A2)
p_complete = formatted.copy()

# Merge with GBMI for #CHR and POS
p_complete = p_complete.merge(
    eur[["rsid", "#CHR", "POS"]], left_on="SNP", right_on="rsid", how="left"
)
p_complete.drop(columns=["rsid"], inplace=True)

# Fill remaining #CHR/POS from TAGC
p_complete = p_complete.merge(
    tagc[["rsid", "chr", "position"]], left_on="SNP", right_on="rsid", how="left"
)
p_complete["#CHR"] = p_complete["#CHR"].fillna(p_complete["chr"])
p_complete["POS"] = p_complete["POS"].fillna(p_complete["position"])
p_complete.drop(columns=["rsid", "chr", "position"], inplace=True)

# Compute p-value (two-sided) and -log10(p)
p_complete["p_value"] = 2 * norm.sf(np.abs(p_complete["Z"]))
p_complete["-log10(P-Value)"] = -np.log10(p_complete["p_value"])

# Drop rows missing positions
before = len(p_complete)
p_complete = p_complete.dropna(subset=["#CHR", "POS"])
print(f"  Dropped {before - len(p_complete)} SNPs without positions")

# Select and order final columns
p_complete = p_complete[
    ["SNP", "A1", "A2", "Z", "#CHR", "POS", "BETA", "SE", "p_value", "-log10(P-Value)"]
]

p_complete.to_csv(OUTPUT_P_COMPLETE, sep="\t", index=False, compression="gzip")
print(f"\nSaved {OUTPUT_P_COMPLETE}: {len(p_complete)} SNPs")
print("Done.")
