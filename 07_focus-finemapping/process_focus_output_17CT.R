# Load the necessary libraries
library(data.table)
# Load necessary library
library(glmnet)
#library(plink2R)
# Load ggplot2
library(ggplot2)
library(dplyr)



focus_file <- "${PROJECT_DIR}/FOCUS/focus_credset_gene_tissue_17CT.tsv"

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


# No PIP threshold for the cell-type panel: we retain the full 90% credible set.


# Remove duplicate gene–tissue pairs
focus_unique <- focus_dt %>%
  distinct(CHR, gene, tissues)


# Map Ensembl IDs (from FOCUS) to gene symbols (used by 1k1k weights/transcripts)
annot <- fread("${PROJECT_DIR}/FOCUS/gene_annotation.txt.gz",
               header = FALSE, sep = "\t")
annot_map <- unique(annot[, .(ensg = sub("\\..*", "", V7), symbol = V4)])

n_before <- nrow(focus_unique)
focus_unique <- focus_unique %>%
  left_join(annot_map, by = c("gene" = "ensg")) %>%
  filter(!is.na(symbol)) %>%
  mutate(gene = symbol) %>%
  select(-symbol)
cat("ENSG->symbol mapping: kept", nrow(focus_unique), "of", n_before, "rows\n")


twas_path <- "${PROJECT_DIR}/TWAS/rebuilt_ct/twas_path_table_ct_new.csv"

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


out_dir <- "${PROJECT_DIR}/gene_model/score_function/focus_geneList_17CT/"
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
