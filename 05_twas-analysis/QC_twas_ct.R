library(data.table)

# ============================================================
# QC for cell-type TWAS results
#   Keep gene-CT pairs with:
#     1. HSQ > 0
#     2. Nonzero predictor variance (TWAS.Z not NA)
#     3. MODELCV.R2 > 0
#
# Writes:
#   - twas_qc_genelist_ct.csv    : one row per passing gene-CT pair
#                                   columns: PANEL, gene_id
#   - twas_qc_summary_ct.csv     : per-CT counts (total / HSQ>0 / nzVar / R2>0 / all3)
# ============================================================

twas_dir     <- "${PROJECT_DIR}/TWAS/TWAS_results_ct"
pos_dir      <- "${PROJECT_DIR}/TWAS/TWAS_pos_ct"
genelist_out <- "${PROJECT_DIR}/TWAS/twas_qc_genelist_ct.csv"
summary_out  <- "${PROJECT_DIR}/TWAS/twas_qc_summary_ct.csv"

# Cell types from pos files
pos_files <- list.files(pos_dir, pattern = "^twas_.*\\.pos$", full.names = TRUE)
tissues   <- gsub("^twas_|\\.pos$", "", basename(pos_files))
cat("Found", length(tissues), "cell types\n\n")

pass_rows    <- list()
summary_rows <- list()

for (tissue in tissues) {
  all_file <- file.path(twas_dir, paste0("TWAS_metaGWAS_", tissue, "_all.dat"))
  if (file.exists(all_file)) {
    dat <- fread(all_file, header = TRUE)
  } else {
    chr_files <- list.files(twas_dir,
                            pattern = paste0("^TWAS_metaGWAS_", tissue, "_chr[0-9]+\\.dat$"),
                            full.names = TRUE)
    if (length(chr_files) == 0) {
      message("WARNING: no TWAS files for cell type: ", tissue); next
    }
    dat <- rbindlist(lapply(chr_files, fread, header = TRUE))
  }

  n_total <- nrow(dat)
  n_hsq   <- sum(!is.na(dat$HSQ) & dat$HSQ > 0)
  n_var   <- sum(!is.na(dat$TWAS.Z))
  n_r2    <- sum(!is.na(dat$MODELCV.R2) & dat$MODELCV.R2 > 0)

  pass <- dat[!is.na(HSQ) & HSQ > 0 &
              !is.na(TWAS.Z) &
              !is.na(MODELCV.R2) & MODELCV.R2 > 0]

  n_pass <- nrow(pass)
  cat(sprintf("%-20s  total=%5d  HSQ>0=%5d  nzVar=%5d  R2>0=%5d  ALL3=%5d\n",
              tissue, n_total, n_hsq, n_var, n_r2, n_pass))

  if (n_pass > 0) {
    pass_rows[[length(pass_rows) + 1]] <- data.table(
      PANEL   = tissue,
      gene_id = pass$ID
    )
  }
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    PANEL = tissue,
    total = n_total, HSQ_gt0 = n_hsq, nzVar = n_var, R2_gt0 = n_r2, ALL3 = n_pass
  )
}

genelist <- rbindlist(pass_rows)
summary  <- rbindlist(summary_rows)

fwrite(genelist, genelist_out)
fwrite(summary,  summary_out)

cat("\n========================================\n")
cat("Total passing gene-CT pairs:", nrow(genelist), "\n")
cat("Unique gene symbols:        ", length(unique(genelist$gene_id)), "\n")
cat("Wrote:", genelist_out, "\n")
cat("Wrote:", summary_out,  "\n")
