# Load the necessary libraries
library(data.table)
# Load necessary library
library(glmnet)
#library(plink2R)
# Load ggplot2
library(ggplot2)
library(dplyr)



focus_file <- "${PROJECT_DIR}/FOCUS/focus_credset_pip_gene_tissue_39.tsv"

focus_dt <- fread(focus_file)

# Quick sanity check
dim(focus_dt)
head(focus_dt)


# Create consistent names
focus_dt <- focus_dt %>%
  rename(
    CHR   = chrom,
    gene  = symbol,
  )


# keep only credible-set genes
focus_dt <- focus_dt %>%
  filter(in_cred_set_pop1 == 1)


# Remove duplicate gene–tissue pairs
focus_unique <- focus_dt %>%
  distinct(CHR, gene, tissues)



twas_path <- "${PROJECT_DIR}/TWAS/rebuilt/twas_path_table_39_new.csv"

twas_path_table <- fread(twas_path)

# Sanity check
#head(twas_path_table)
#colnames(twas_path_table)


focus_with_paths <- focus_unique %>%
  left_join(
    twas_path_table,
    by = c("tissues" = "PANEL")
  )

focus_final <- focus_with_paths %>%
  select(CHR, gene, tissue)

print(anyNA(focus_final))


out_dir <- "${PROJECT_DIR}/gene_model/score_function/focus_geneList_39/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Unique chromosome–tissue combinations
unique_combos <- focus_final %>%
  distinct(CHR, tissue)

for (i in seq_len(nrow(unique_combos))) {
  chr <- unique_combos$CHR[i]
  t   <- unique_combos$tissue[i]

  subset_df <- focus_final %>%
    filter(CHR == chr, tissue == t)

  filename <- paste0(
    "fdr_sig_genes_chr", chr, "_",
    gsub("[ /]", "_", t),
    ".txt"
  )

  writeLines(
    unique(subset_df$gene),
    file.path(out_dir, filename)
  )
}
