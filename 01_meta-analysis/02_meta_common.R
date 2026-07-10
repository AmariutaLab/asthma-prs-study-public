# ============================================================
# Stage 2: Fixed-effects meta-analysis of common SNPs
# Input:  TAGC EUR + GBMI EUR raw summary statistics
# Output: output/meta_analysis_results_common.csv
#         (inverse-variance weighted combined beta/SE for
#          SNPs present in both TAGC and GBMI)
# ============================================================

library(data.table)
library(dplyr)
library(metafor)
library(tidyr)
library(purrr)
library(progress)

# --- Configuration ---
tagc_path <- "../TAGC/TAGC meta-analyses results for asthma risk/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv"
gbmi_path <- "../GBMI/sumstats/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz"
output_path <- "output/meta_analysis_results_common.csv"

# --- Read data ---
tagc <- fread(tagc_path, sep = "\t")
eur <- fread(cmd = paste("gunzip -c", gbmi_path), sep = "\t")

# --- Per-study QC: drop variants with missing or incomplete effect-size / SE ---
n_tagc_before <- nrow(tagc)
tagc <- tagc[!is.na(European_ancestry_beta_fix) &
             !is.na(European_ancestry_se_fix) &
             European_ancestry_se_fix > 0]
cat("TAGC: dropped", n_tagc_before - nrow(tagc),
    "variants with missing/zero beta or SE (", n_tagc_before, "->", nrow(tagc), ")\n")

n_eur_before <- nrow(eur)
eur <- eur[!is.na(inv_var_meta_beta) &
           !is.na(inv_var_meta_sebeta) &
           inv_var_meta_sebeta > 0]
cat("GBMI: dropped", n_eur_before - nrow(eur),
    "variants with missing/zero beta or SE (", n_eur_before, "->", nrow(eur), ")\n")

# --- Inner join on common SNPs ---
snps <- tagc %>%
  select(chr, rsid, reference_allele, alternate_allele,
         European_ancestry_beta_fix, European_ancestry_se_fix) %>%
  inner_join(
    eur %>% select(`#CHR`, rsid, REF, ALT, inv_var_meta_beta, inv_var_meta_sebeta),
    by = c("chr" = "#CHR", "rsid" = "rsid")
  )

cat("Common SNPs found:", nrow(snps), "\n")

# --- Allele concordance check and harmonization ---
# TAGC beta is relative to alternate_allele; GBMI beta is relative to ALT.
# Check whether effect alleles agree, and flip GBMI beta if reversed.
snps <- snps %>%
  mutate(
    allele_status = case_when(
      alternate_allele == ALT & reference_allele == REF ~ "concordant",
      alternate_allele == REF & reference_allele == ALT ~ "flipped",
      TRUE ~ "mismatch"
    )
  )

n_concordant <- sum(snps$allele_status == "concordant")
n_flipped <- sum(snps$allele_status == "flipped")
n_mismatch <- sum(snps$allele_status == "mismatch")
cat("Allele concordance:\n")
cat("  Concordant:", n_concordant, "\n")
cat("  Flipped:   ", n_flipped, "(GBMI beta sign will be reversed)\n")
cat("  Mismatch:  ", n_mismatch, "(removed — alleles incompatible)\n")

# Flip GBMI beta for reversed alleles
snps <- snps %>%
  mutate(
    inv_var_meta_beta = ifelse(allele_status == "flipped",
                               -inv_var_meta_beta,
                               inv_var_meta_beta)
  ) %>%
  filter(allele_status != "mismatch") %>%
  select(-reference_allele, -alternate_allele, -REF, -ALT, -allele_status) %>%
  rename(
    beta_tagc = European_ancestry_beta_fix,
    se_tagc = European_ancestry_se_fix,
    beta_gbmi = inv_var_meta_beta,
    se_gbmi = inv_var_meta_sebeta
  )

cat("SNPs after harmonization:", nrow(snps), "\n")

# --- Fixed-effects inverse-variance weighted meta-analysis ---
perform_meta_analysis <- function(beta_tagc, se_tagc, beta_gbmi, se_gbmi,
                                  progress_bar = NULL) {
  res <- rma(yi = c(beta_tagc, beta_gbmi),
             sei = c(se_tagc, se_gbmi),
             method = "FE")
  if (!is.null(progress_bar)) progress_bar$tick()
  return(c(res$beta, res$se))
}

pb <- progress_bar$new(
  format = "  Progress [:bar] :percent in :elapsed (ETA: :eta)",
  total = nrow(snps), clear = FALSE
)

meta_results <- snps %>%
  mutate(
    meta_analysis = pmap(
      list(beta_tagc, se_tagc, beta_gbmi, se_gbmi),
      ~ perform_meta_analysis(..1, ..2, ..3, ..4, progress_bar = pb)
    )
  ) %>%
  mutate(
    combined_beta = map_dbl(meta_analysis, 1),
    combined_se = map_dbl(meta_analysis, 2)
  ) %>%
  select(-meta_analysis)

# --- Save ---
write.csv(meta_results, output_path, row.names = FALSE)
cat("Saved meta-analysis results to:", output_path, "\n")
