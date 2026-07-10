# ============================================================
# Build the gene-tissue table consumed downstream by FOCUS (39 GTEx tissues).
# For every tissue in TissueGroups.txt this script:
#   1. Loads the per-tissue TWAS .dat files from TWAS_results/
#   2. Restricts to genes present in BOTH TranscriptsIn<Tissue>Model.txt and
#      TranscriptsIn<Tissue>Model_keep.txt (TCSC heritable-genes filter)
#   3. Keeps only rows with TWAS.P < 0.05
#   4. Joins gene-annotation columns (CHR, START, STOP, SYM) and writes the
#      combined, p-value-sorted table to FOCUS/sorted_gene_all39Tissue_new.csv.
#
# Output schema matches the legacy sorted_gene format (gene_id, V1=TWAS.Z,
# p_value, PANEL, weight, WGT, N, ID, CHR, P0, P1, SYM, START, STOP) so it can
# drop into the existing FOCUS pipeline.
# ============================================================

library(data.table)

# ============================================================
# Paths — adjust as needed
# ============================================================
twas_dir       <- "${PROJECT_DIR}/TWAS/TWAS_results"
tissue_file    <- "${REF_DIR}/TCSC/analysis/TissueGroups.txt"
gene_anno_file <- "${REF_DIR}/gene_annotation.txt.gz"
herit_dir      <- "${REF_DIR}/TCSC/weights/heritablegenes"
output_file    <- "${PROJECT_DIR}/FOCUS/sorted_gene_all39Tissue_new.csv"

# ============================================================
# 1. Read tissue groups (same file used in TCSC pipeline)
# ============================================================
y <- fread(tissue_file, header = TRUE)
tissues    <- unique(y$MetaTissue)
n_eqtl     <- sapply(tissues, function(t) sum(y$N_EUR[y$MetaTissue == t]))
small_flag <- n_eqtl < 320   # TRUE = small tissue (Nall); FALSE = subsampled to 320

cat("Found", length(tissues), "unique MetaTissues\n")

# ============================================================
# 2. Read gene annotations
# ============================================================
genes_full <- fread(gene_anno_file, sep = "\t", header = FALSE)
genes <- genes_full[, .(V1, V2, V3, V4, V7)]
setnames(genes, c("CHR", "START", "STOP", "ID", "SYM"))
# Remove version from ENSG ID
genes[, SYM := sub("\\..*", "", SYM)]

cat("Gene annotation loaded:", nrow(genes), "entries\n")

# ============================================================
# 3. For each tissue, read TWAS results and filter
# ============================================================
all_results <- list()
total_loaded <- 0
total_heritable <- 0
total_sig <- 0

for (i in seq_along(tissues)) {
  tissue       <- tissues[i]
  tissue_lower <- tolower(tissue)

  # Try _all.dat first
  all_file <- file.path(twas_dir, paste0("twas_", tissue_lower, "_all.dat"))

  if (file.exists(all_file)) {
    dat <- fread(all_file, header = TRUE)
  } else {
    # Fall back: combine per-chromosome files
    chr_files <- list.files(twas_dir,
                            pattern = paste0("^twas_", tissue_lower, "_chr[0-9]+\\.dat$"),
                            full.names = TRUE)
    if (length(chr_files) == 0) {
      message("WARNING: No TWAS files found for tissue: ", tissue, " (expected pattern: twas_", tissue_lower, "_chr*.dat)")
      next
    }
    dat <- rbindlist(lapply(chr_files, fread, header = TRUE))
  }

  total_loaded <- total_loaded + nrow(dat)
  cat(tissue, ":", nrow(dat), "genes loaded,")

  # Extract ENSG ID (with version) from the FILE column
  # FILE example: .../Heart_Atrial_Appendage.ENSG00000228794.8.wgt.RDat
  dat[, gene_id := sub("^.*\\.(ENSG[0-9]+\\.[0-9]+)\\.wgt\\.RDat$", "\\1", FILE)]

  # ------------------------------------------------------------------
  # Filter to genes present in BOTH TranscriptsIn{Tissue}Model.txt
  # AND TranscriptsIn{Tissue}Model_keep.txt
  # (same logic as Make_Pos_File_forWeights_multi.R)
  # ------------------------------------------------------------------
  prefix <- ifelse(small_flag[i], "Nall", "N320")
  transcripts_file <- file.path(herit_dir, prefix, paste0("TranscriptsIn", tissue, "Model.txt"))
  keep_file        <- file.path(herit_dir, prefix, paste0("TranscriptsIn", tissue, "Model_keep.txt"))

  if (file.exists(transcripts_file) && file.exists(keep_file)) {
    transcripts_all  <- fread(transcripts_file, header = FALSE)$V1
    transcripts_keep <- fread(keep_file, header = FALSE)$V1
    valid_genes <- intersect(transcripts_all, transcripts_keep)
    dat <- dat[gene_id %in% valid_genes]
    total_heritable <- total_heritable + nrow(dat)
    cat(nrow(dat), "in both transcript files,")
  } else {
    if (!file.exists(transcripts_file)) message("WARNING: transcripts file not found: ", transcripts_file)
    if (!file.exists(keep_file))        message("WARNING: keep file not found: ", keep_file)
    message("  — skipping heritability filter for ", tissue)
  }

  # Filter by TWAS.P < 0.05
  dat <- dat[TWAS.P < 0.05, ]
  total_sig <- total_sig + nrow(dat)
  cat(nrow(dat), "pass p < 0.05\n")

  if (nrow(dat) == 0) next

  # Build columns matching the original sorted_gene format
  dat[, V1      := TWAS.Z]
  dat[, p_value := TWAS.P]
  dat[, PANEL   := tissue_lower]
  dat[, weight  := FILE]
  dat[, WGT     := sub("${REF_DIR}/gtex/weights/", "", FILE)]
  dat[, N       := ifelse(small_flag[i], n_eqtl[i], 320)]

  all_results[[length(all_results) + 1]] <- dat[, .(gene_id, V1, p_value, PANEL, weight, WGT, N)]
}

cat("\n===== TOTAL SUMMARY =====\n")
cat("Total genes loaded:              ", total_loaded, "\n")
cat("Total in both transcript files:  ", total_heritable, "\n")
cat("Total pass p < 0.05:             ", total_sig, "\n")
cat("=========================\n\n")

# ============================================================
# 4. Combine all tissues and sort by p-value
# ============================================================
sorted_gene <- rbindlist(all_results, fill = TRUE)
sorted_gene <- sorted_gene[order(p_value)]
cat("\nTotal gene-tissue pairs (p < 0.05):", nrow(sorted_gene), "\n")

# ============================================================
# 5. Add gene annotation columns
# ============================================================
# Match by ENSG ID (remove version first)
transcripts <- sub("\\..*", "", sorted_gene$gene_id)
m <- match(transcripts, genes$SYM)

sorted_gene[, ID    := genes$ID[m]]        # gene symbol (e.g., AGRN)
sorted_gene[, CHR   := genes$CHR[m]]       # e.g., chr1
sorted_gene[, P0    := genes$START[m] - 1000000]  # buffer for genome build
sorted_gene[, P1    := genes$STOP[m]  + 1000000]
sorted_gene[, SYM   := genes$SYM[m]]       # ENSG ID without version
sorted_gene[, START := genes$START[m]]
sorted_gene[, STOP  := genes$STOP[m]]

# Remove version from ID (gene symbol — usually no-op, but matches original script)
sorted_gene[, ID := sub("\\..*", "", ID)]

# Append PANEL to SYM (matches original format: ENSG00000228794_heart_atrial_appendage)
sorted_gene[, SYM := paste0(SYM, "_", PANEL)]

# Report unmatched genes
n_na <- sum(is.na(sorted_gene$CHR))
if (n_na > 0) {
  cat("WARNING:", n_na, "gene-tissue pairs could not be matched to gene_annotation.txt.gz\n")
  cat("  These will have NA in CHR/P0/P1/START/STOP columns\n")
}

# ============================================================
# 6. Write output
# ============================================================
write.table(sorted_gene, file = output_file,
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
cat("Written to:", output_file, "\n")

