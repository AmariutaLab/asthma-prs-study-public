log_file <- "ld_clumping.log"
sink(log_file, append = TRUE, split = TRUE)


# Load the necessary libraries
library(data.table)
library(glmnet)


# ==========================================================================
# Paths
# ==========================================================================
path_table_file <- "${PROJECT_DIR}/TWAS/rebuilt/twas_path_table_39_new.csv"
tissue_file     <- "${REF_DIR}/TCSC/analysis/TissueGroups.txt"
gene_anno_file  <- "${REF_DIR}/gene_annotation.txt.gz"

# Predicted-expression matrices (unchanged — same GTEx v8 weights as before)
pred_base       <- "${REF_DIR}/TCSC/predicted_expression"

# Original transcripts files (these correspond to designmat columns; the
# rebuilt transcripts files are row-aligned with the NEW gtex_twas Z-scores
# but NOT with designmat, so we need both)
orig_trans_base <- "${REF_DIR}/TCSC/weights/heritablegenes"

# Output
output_dir      <- "twas_gene_lists"


# ==========================================================================
# Load gene annotations
# ==========================================================================
genes_full <- fread(gene_anno_file, sep = "\t", header = FALSE)
genes <- genes_full[, .(V1, V2, V3, V4, V7)]
setnames(genes, old = c("V1", "V2", "V3", "V4", "V7"),
                new = c("CHR", "START", "STOP", "ID", "SYM"))

# remove version
genes$SYM <- sapply(seq_along(genes$SYM),
                   function(x) strsplit(genes$SYM[x], split = "[.]")[[1]][1])
genes[, CHR := as.character(CHR)]


# ==========================================================================
# Tissue groups (only used to decide N320 vs Nall prefix for pred + orig
# transcripts paths; the rebuilt path table already encodes the new TWAS
# and keep-file locations)
# ==========================================================================
y <- fread(tissue_file, header = TRUE)
tissues_grp <- unique(y$MetaTissue)
n_eqtl      <- sapply(seq_along(tissues_grp),
                      function(x) sum(y$N_EUR[which(y$MetaTissue == tissues_grp[x])]))
small_set   <- tissues_grp[n_eqtl < 320]


# ==========================================================================
# Build file_sets from the rebuilt path table
# ==========================================================================
path_table <- fread(path_table_file)

file_sets <- list()
for (i in seq_len(nrow(path_table))) {
  row <- path_table[i]
  tissue <- row$tissue

  if (tissue %in% small_set) {
    prefix    <- "Nall"
    pred_mid  <- paste0(tissue, "_v8_allEUR")
  } else {
    prefix    <- "N320"
    pred_mid  <- paste0(tissue, "_v8_320EUR")
  }

  pred_file       <- file.path(pred_base, prefix,
                               paste0("designmat_", pred_mid, "_double.RData"))
  orig_trans_file <- file.path(orig_trans_base, prefix,
                               paste0("TranscriptsIn", tissue, "Model.txt"))

  file_sets[[length(file_sets) + 1]] <- list(
    transcripts          = row$transcripts,           # NEW: row-aligned with gtex_twas
    gtex_twas            = row$gtex_twas,             # NEW: Marginal_alphas_NEW_TWAS_*.txt.gz
    transcript_keep_file = row$transcript_keep_file,  # original heritability keep list
    pred                 = pred_file,                 # original GTEx designmat
    transcripts_pred     = orig_trans_file,           # original transcripts (matches designmat)
    pos_df_col1          = tissue
  )
}


# ==========================================================================
# LD clumping
# ==========================================================================
ld_clump_genes <- function(sorted_gene, corr_mat, r2_threshold = 0.1) {

  # sorted_gene must already be sorted by p_value ascending
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


process_tissue <- function(tissue_entry, p_thresholds) {

  tissue_name <- tissue_entry$pos_df_col1
  cat("\n==============================\n")
  cat("Processing tissue:", tissue_name, "\n")
  cat("==============================\n")

  # --------------------------------------
  # 1. Predicted expression -> corr matrix
  # --------------------------------------
  load(tissue_entry$pred)   # loads df1
  corr_matrix <- cor(df1)

  # Label corr_matrix using the ORIGINAL transcripts file (row-aligned with df1)
  transcripts_pred <- fread(tissue_entry$transcripts_pred, header = FALSE)
  setnames(transcripts_pred, "V1", "gene_id")
  corr_base <- sub("\\..*$", "", transcripts_pred$gene_id)
  colnames(corr_matrix) <- corr_base
  rownames(corr_matrix) <- corr_base

  # --------------------------------------
  # 2. NEW transcripts (row-aligned with new gtex_twas)
  # --------------------------------------
  transcripts <- fread(tissue_entry$transcripts, header = FALSE)
  setnames(transcripts, "V1", "gene_id")
  transcripts$gene_base <- sub("\\..*$", "", transcripts$gene_id)

  # --------------------------------------
  # 3. Heritability keep-list (unchanged)
  # --------------------------------------
  transcripts_keep <- fread(tissue_entry$transcript_keep_file, header = FALSE)
  transcripts_keep$gene_base <- sub("\\..*$", "", transcripts_keep$V1)

  # --------------------------------------
  # 4. NEW TWAS Z-scores
  # --------------------------------------
  gtex_twas <- fread(tissue_entry$gtex_twas, header = FALSE)
  gtex_combine <- cbind(transcripts, gtex_twas)
  gtex_combine[, p_value := 2 * pnorm(abs(V1), lower.tail = FALSE)]

  # rank by p-value ascending
  sorted_gene <- gtex_combine[order(gtex_combine$p_value), ]
  sorted_gene$gene_base <- sub("\\..*$", "", sorted_gene$gene_id)
  sorted_gene$gene_id   <- sorted_gene$gene_base

  sorted_gene$gene_id <- trimws(as.character(sorted_gene$gene_id))
  rownames(corr_matrix) <- trimws(rownames(corr_matrix))
  colnames(corr_matrix) <- trimws(colnames(corr_matrix))

  # --------------------------------------
  # 5. LD Clumping
  # --------------------------------------
  clumped_genes <- ld_clump_genes(sorted_gene, corr_matrix)

  # Filter through transcripts_keep
  clumped_filtered <- clumped_genes[clumped_genes %in% transcripts_keep$gene_base]

  cat("Clumped genes kept:", length(clumped_filtered), "\n")

  # --------------------------------------
  # 6. Thresholding and saving lists
  # --------------------------------------
  for (thr in p_thresholds) {
    df_thr <- subset(
      sorted_gene,
      gene_id %in% clumped_filtered & p_value < thr
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
    file_name <- paste0(output_dir, "/", tissue_name, "_gene_list_p", thr_label, ".txt")

    write.table(
      df_map,
      file_name,
      quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    cat("Saved:", file_name, "(", length(gene_list), "genes )\n")
  }
  cat("Finished:", tissue_name, "\n")
}


# ==========================================================================
# Run
# ==========================================================================
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p_thresholds <- c(5e-5, 5e-4, 5e-3, 0.05)

for (i in seq_along(file_sets)) {
  process_tissue(file_sets[[i]], p_thresholds)
}


sink()
