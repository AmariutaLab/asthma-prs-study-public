log_file <- "ld_clumping.log"
sink(log_file, append = TRUE, split = TRUE)


# Load the necessary libraries
library(data.table)


# ==========================================================================
# Paths
# ==========================================================================
path_table_file <- "${PROJECT_DIR}/TWAS/rebuilt_ct/twas_path_table_ct_new.csv"
gene_anno_file  <- "${REF_DIR}/gene_annotation.txt.gz"

# Cell-type predicted-expression matrices (collaborator)
pred_base       <- "${REF_DIR}/1k1k_cell_type/celltype_pred_expr"

# Output
output_dir      <- "twas_gene_lists"


# ==========================================================================
# Load gene annotations
# Columns: V1=chr (with "chr" prefix), V2=start, V3=stop, V4=symbol, V7=Ensembl
# CT gene IDs are symbols → merge on V4.
# ==========================================================================
genes_full <- fread(gene_anno_file, sep = "\t", header = FALSE)
genes <- genes_full[, .(V1, V2, V3, V4, V7)]
setnames(genes, c("CHR", "START", "STOP", "SYM", "ENSG"))
genes[, CHR := as.character(CHR)]


# ==========================================================================
# Build entries from CT path table
# CT path table fields: tissue, gtex_twas, transcripts, weights, PANEL
# (no transcript_keep_file; QC already enforced during rebuild)
# weights path basename gives case-preserved CT name used by collaborator CSVs.
# ==========================================================================
path_table <- fread(path_table_file)

file_sets <- list()
for (i in seq_len(nrow(path_table))) {
  row <- path_table[i]

  ct_dir <- basename(sub("/$", "", row$weights))
  pred_file <- file.path(pred_base, paste0(ct_dir, "_pred_expr.csv"))

  file_sets[[length(file_sets) + 1]] <- list(
    transcripts = row$transcripts,
    gtex_twas   = row$gtex_twas,
    pred        = pred_file,
    pos_name    = ct_dir
  )
}


# ==========================================================================
# LD clumping
# ==========================================================================
ld_clump_genes <- function(sorted_gene, corr_mat, r2_threshold = 0.1) {

  gene_order <- sorted_gene$gene_id
  selected <- c()

  for (g in gene_order) {
    if (!(g %in% rownames(corr_mat))) {
      cat("\n Missing ROW in corr_mat:", g, "\n")
      next
    }
    if (!(g %in% colnames(corr_mat))) {
      cat("\n Missing COL in corr_mat:", g, "\n")
      next
    }

    if (length(selected) == 0) {
      selected <- c(selected, g)
      next
    }

    r_vec  <- corr_mat[g, selected]
    r2_vec <- r_vec^2

    if (all(r2_vec < r2_threshold, na.rm = TRUE)) {
      selected <- c(selected, g)
    }
  }

  return(selected)
}


process_ct <- function(entry, p_thresholds) {

  ct_name <- entry$pos_name
  cat("\n==============================\n")
  cat("Processing CT:", ct_name, "\n")
  cat("==============================\n")

  # --------------------------------------
  # 1. Predicted expression -> corr matrix
  # CSV layout: row 1 = header (empty cell + gene symbols), col 1 = sample ID
  # --------------------------------------
  pe <- fread(entry$pred)
  pe[, 1 := NULL]                      # drop sample-ID column
  mat <- as.matrix(pe)
  corr_matrix <- cor(mat)
  # corr_matrix already has gene symbols as row/col names from data.table

  # --------------------------------------
  # 2. NEW transcripts (row-aligned with new gtex_twas)
  # --------------------------------------
  transcripts <- fread(entry$transcripts, header = FALSE)
  setnames(transcripts, "V1", "gene_id")

  # --------------------------------------
  # 3. NEW TWAS Z-scores
  # --------------------------------------
  gtex_twas <- fread(entry$gtex_twas, header = FALSE)
  gtex_combine <- cbind(transcripts, gtex_twas)
  gtex_combine[, p_value := 2 * pnorm(abs(V1), lower.tail = FALSE)]

  sorted_gene <- gtex_combine[order(p_value), ]
  sorted_gene$gene_id <- trimws(as.character(sorted_gene$gene_id))
  rownames(corr_matrix) <- trimws(rownames(corr_matrix))
  colnames(corr_matrix) <- trimws(colnames(corr_matrix))

  # --------------------------------------
  # 4. LD Clumping
  # --------------------------------------
  clumped_genes <- ld_clump_genes(sorted_gene, corr_matrix)
  cat("Clumped genes kept:", length(clumped_genes), "\n")

  # --------------------------------------
  # 5. Thresholding and saving lists
  # --------------------------------------
  for (thr in p_thresholds) {
    df_thr <- subset(
      sorted_gene,
      gene_id %in% clumped_genes & p_value < thr
    )

    gene_list <- df_thr$gene_id

    df_map <- merge(
      data.frame(gene_id = gene_list),
      genes,
      by.x = "gene_id",
      by.y = "SYM",
      all.x = TRUE
    )

    thr_label <- gsub("\\.", "_", as.character(thr))
    file_name <- paste0(output_dir, "/", ct_name, "_gene_list_p", thr_label, ".txt")

    write.table(
      df_map,
      file_name,
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    cat("Saved:", file_name, "(", length(gene_list), "genes )\n")
  }
  cat("Finished:", ct_name, "\n")
}


# ==========================================================================
# Run
# ==========================================================================
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p_thresholds <- c(5e-5, 5e-4, 5e-3, 0.05)

for (i in seq_along(file_sets)) {
  process_ct(file_sets[[i]], p_thresholds)
}


sink()
