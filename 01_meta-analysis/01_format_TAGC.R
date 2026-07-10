# ============================================================
# Stage 1: Format TAGC summary statistics
# Input:  Raw TAGC sumstats (.tsv)
# Output: Standardized columns (A1=ALT, A2=REF, SNP, N, Z)
# ============================================================

library(dplyr)
library(data.table)

# --- Configuration ---
ancestry <- "European_ancestry"  # Options: European_ancestry, Multiancestry
method <- "fix"                  # Options: fix (fixed-effects), rand (random-effects)
input_path <- "../TAGC/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv"
output_path <- paste0("output/TAGC_formatted_", ancestry, "_", method, ".txt.gz")

# --- Read data ---
tagc <- fread(input_path, sep = "\t")

# --- Validate columns ---
beta_col <- paste0(ancestry, "_beta_", method)
se_col <- paste0(ancestry, "_se_", method)

if (!all(c(beta_col, se_col) %in% colnames(tagc))) {
  stop("Columns not found: ", beta_col, ", ", se_col,
       ". Check ancestry and method settings.")
}

# Sample sizes (from TAGC publication)
N_value <- ifelse(ancestry == "European_ancestry", 127669,
                  ifelse(ancestry == "Multiancestry", 142486, NA))
if (is.na(N_value)) stop("Invalid ancestry.")

# --- QC: drop variants with missing or incomplete effect-size / standard-error info ---
n_before <- nrow(tagc)
tagc <- tagc[!is.na(tagc[[beta_col]]) & !is.na(tagc[[se_col]]) & tagc[[se_col]] > 0, ]
cat("Dropped", n_before - nrow(tagc),
    "TAGC variants with missing/zero beta or SE (", n_before, "->", nrow(tagc), ")\n")

# --- Transform ---
# A1 = alternate_allele = effect allele (ALT)
# A2 = reference_allele = reference allele (REF)
# Beta in TAGC is relative to alternate_allele, so Z = beta/SE is relative to A1
tagc_formatted <- tagc %>%
  mutate(
    A1 = alternate_allele,
    A2 = reference_allele,
    SNP = rsid,
    N = N_value,
    Z = get(beta_col) / get(se_col)
  ) %>%
  select(A1, A2, SNP, N, Z)

# --- Save ---
fwrite(tagc_formatted, file = output_path, compress = "gzip", sep = "\t")
cat("Saved formatted TAGC to:", output_path, "\n")
cat("Total SNPs:", nrow(tagc_formatted), "\n")
