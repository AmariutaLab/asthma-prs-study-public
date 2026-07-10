library(data.table)

# ============================================================
# Build .pos files for OneK1K cell-type weights (v2)
#   Weights:     ${REF_DIR}/1k1k_cell_type/celltype_weights_v2/{ct}/{gene}.wgt.RDat
#   Annotation:  ${PROJECT_DIR}/FOCUS/gene_annotation.txt.gz
#   N per CT:    number of donors in celltype_eQTL_v2/{ct}/plinkData/*.fam
# ============================================================

weights_base <- "${REF_DIR}/1k1k_cell_type/celltype_weights_v2"
eqtl_base    <- "${REF_DIR}/1k1k_cell_type/celltype_eQTL_v2"
annot_file   <- "${PROJECT_DIR}/FOCUS/gene_annotation.txt.gz"
outdir       <- "${PROJECT_DIR}/TWAS/TWAS_pos_ct"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Gene annotation: V1=chr, V2=start, V3=stop, V4=symbol, V7=ENSG (versioned)
# ------------------------------------------------------------
genes_full <- fread(annot_file, sep = "\t", header = FALSE)
genes <- genes_full[, .(V1, V2, V3, V4, V7)]
setnames(genes, c("CHR", "START", "STOP", "SYM", "ENSG"))
genes[, ENSG := sub("\\..*", "", ENSG)]
genes <- unique(genes, by = "SYM")  # in case of duplicate symbols, keep first

# ------------------------------------------------------------
# Iterate cell types
# ------------------------------------------------------------
cts <- list.dirs(weights_base, recursive = FALSE, full.names = FALSE)
cat("Found", length(cts), "cell types\n")

for (ct in cts) {
  wgt_files <- list.files(file.path(weights_base, ct),
                          pattern = "\\.wgt\\.RDat$", full.names = FALSE)
  if (length(wgt_files) == 0) {
    message("WARNING: no weight files for ", ct); next
  }

  # N = donors in any fam for this cell type (same across genes within a ct)
  fam_files <- list.files(file.path(eqtl_base, ct, "plinkData"),
                          pattern = "\\.fam$", full.names = TRUE)
  if (length(fam_files) == 0) {
    message("WARNING: no fam files for ", ct, " — setting N = NA"); n_ct <- NA_integer_
  } else {
    n_ct <- nrow(fread(fam_files[1], header = FALSE))
  }

  gene_syms <- sub("\\.wgt\\.RDat$", "", wgt_files)
  dt <- data.table(SYM = gene_syms,
                   WGT = paste0(ct, "/", wgt_files))

  m <- match(dt$SYM, genes$SYM)
  dt[, CHR  := genes$CHR[m]]
  dt[, P0   := genes$START[m] - 1000000]   # 1Mb buffer (annotation may differ from plink build)
  dt[, P1   := genes$STOP[m]  + 1000000]
  dt[, ID   := genes$SYM[m]]                # keep gene symbol as ID (matches .wgt filename)
  dt[, PANEL := tolower(ct)]
  dt[, N    := n_ct]

  # numeric chr only
  dt[, CHR := suppressWarnings(as.numeric(gsub("chr", "", CHR)))]

  missing <- dt[is.na(CHR)]
  if (nrow(missing) > 0) {
    cat(sprintf("  %s: %d / %d genes unmatched in annotation (dropped)\n",
                ct, nrow(missing), nrow(dt)))
  }
  dt <- dt[!is.na(CHR) & CHR > 0]

  dt_out <- dt[, .(PANEL, WGT, ID, CHR, P0, P1, N)]
  setorder(dt_out, CHR, P0)

  outfile <- file.path(outdir, paste0("twas_", tolower(ct), ".pos"))
  fwrite(dt_out, file = outfile, sep = "\t", col.names = TRUE)

  cat(sprintf("%-18s  %5d genes -> %s (N=%s)\n",
              ct, nrow(dt_out), outfile, n_ct))
}

cat("\nDone. .pos files in:", outdir, "\n")
