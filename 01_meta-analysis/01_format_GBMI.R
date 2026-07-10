# ============================================================
# Stage 1: Format GBMI summary statistics
# Input:  Raw GBMI sumstats (.txt.gz)
# Output: Standardized columns (A1=ALT, A2=REF, SNP, N, Z)
# ============================================================

library(dplyr)
library(data.table)

# --- Configuration ---
ancestry <- "eur"  # Options: eur, afr, amr, eas, sas
input_path <- paste0(
  "../GBMI/sumstats/Asthma_Bothsex_", ancestry,
  "_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
)
output_path <- paste0("output/GBMI_formatted_", ancestry, ".txt.gz")

# --- Read data ---
eur <- fread(cmd = paste("gunzip -c", input_path), sep = "\t")

# --- QC: drop variants with missing or incomplete effect-size / standard-error info ---
n_before <- nrow(eur)
eur <- eur[!is.na(inv_var_meta_beta) & !is.na(inv_var_meta_sebeta) & inv_var_meta_sebeta > 0]
cat("Dropped", n_before - nrow(eur),
    "GBMI variants with missing/zero beta or SE (", n_before, "->", nrow(eur), ")\n")

# --- Transform ---
# A1 = ALT = effect allele
# A2 = REF = reference allele
# Beta in GBMI is relative to ALT, so Z = beta/SE is relative to A1
eur_formatted <- eur %>%
  mutate(
    A1 = ALT,
    A2 = REF,
    SNP = rsid,
    N = N_case + N_ctrl,
    Z = inv_var_meta_beta / inv_var_meta_sebeta
  ) %>%
  select(A1, A2, SNP, N, Z)

# --- Save ---
fwrite(eur_formatted, file = output_path, compress = "gzip", sep = "\t")
cat("Saved formatted GBMI to:", output_path, "\n")
cat("Total SNPs:", nrow(eur_formatted), "\n")
