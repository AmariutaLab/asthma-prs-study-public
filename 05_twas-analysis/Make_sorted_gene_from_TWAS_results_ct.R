# ============================================================
# Build the gene-celltype table consumed downstream by FOCUS (17 OneK1K cell
# types). For every CT with a .pos file in TWAS_pos_ct/ this script:
#   1. Loads the per-CT TWAS .dat files from TWAS_results_ct/
#   2. Restricts to gene-CT pairs that passed QC_twas_ct.R
#      (HSQ>0, nonzero predictor variance, MODELCV.R2>0)
#   3. Keeps only rows with TWAS.P < 0.05
#   4. Joins gene-annotation columns (CHR, START, STOP, SYM) and writes the
#      combined, p-value-sorted table to FOCUS/sorted_gene_all17CT_new.csv.
#
# Output schema is the same legacy sorted_gene format as the 39-tissue
# counterpart so it drops into the existing FOCUS pipeline.
# ============================================================

library(data.table)

# ============================================================
# Paths — adjust as needed
# ============================================================
twas_dir       <- "${PROJECT_DIR}/TWAS/TWAS_results_ct"
pos_dir        <- "${PROJECT_DIR}/TWAS/TWAS_pos_ct"
qc_genelist    <- "${PROJECT_DIR}/TWAS/twas_qc_genelist_ct.csv"
gene_anno_file <- "${REF_DIR}/gene_annotation.txt.gz"
output_file    <- "${PROJECT_DIR}/FOCUS/sorted_gene_all17CT_new.csv"

# Load QC gene list (PANEL, gene_id) produced by QC_twas_ct.R
qc <- fread(qc_genelist)
qc[, key := paste0(PANEL, "|", gene_id)]
cat("QC gene list loaded:", nrow(qc), "gene-CT pairs\n")

# ============================================================
# 1. Discover cell types from pos files (replaces TissueGroups.txt)
# ============================================================
pos_files <- list.files(pos_dir, pattern = "^twas_.*\\.pos$", full.names = TRUE)
tissues   <- gsub("^twas_|\\.pos$", "", basename(pos_files))

cat("Found", length(tissues), "cell types\n")

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
# 3. For each cell type, read TWAS results and filter
# ============================================================
all_results <- list()
total_loaded <- 0
total_sig <- 0

for (i in seq_along(tissues)) {
  tissue       <- tissues[i]
  tissue_lower <- tissue   # already lowercase in pos-file naming

  # N comes from the pos file (uniform across genes within a cell type)
  pos    <- fread(file.path(pos_dir, paste0("twas_", tissue_lower, ".pos")), header = TRUE)
  n_ct   <- as.integer(pos$N[1])

  # Try _all.dat first
  all_file <- file.path(twas_dir, paste0("TWAS_metaGWAS_", tissue_lower, "_all.dat"))

  if (file.exists(all_file)) {
    dat <- fread(all_file, header = TRUE)
  } else {
    # Fall back: combine per-chromosome files
    chr_files <- list.files(twas_dir,
                            pattern = paste0("^TWAS_metaGWAS_", tissue_lower, "_chr[0-9]+\\.dat$"),
                            full.names = TRUE)
    if (length(chr_files) == 0) {
      message("WARNING: No TWAS files found for cell type: ", tissue)
      next
    }
    dat <- rbindlist(lapply(chr_files, fread, header = TRUE))
  }

  total_loaded <- total_loaded + nrow(dat)
  cat(tissue, ":", nrow(dat), "genes loaded,")

  # Gene ID = symbol (CT weights are named by symbol; FUSION propagates ID column from pos)
  dat[, gene_id := ID]

  # Keep only gene-CT pairs that passed QC (HSQ>0, nonzero var, MODELCV.R2>0)
  dat <- dat[paste0(tissue_lower, "|", gene_id) %in% qc$key]
  cat(nrow(dat), "pass QC,")

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
  dat[, WGT     := sub("${REF_DIR}/1k1k_cell_type/celltype_weights_v2/", "", FILE)]
  dat[, N       := n_ct]

  all_results[[length(all_results) + 1]] <- dat[, .(gene_id, V1, p_value, PANEL, weight, WGT, N)]
}

cat("\n===== TOTAL SUMMARY =====\n")
cat("Total genes loaded:              ", total_loaded, "\n")
cat("Total pass p < 0.05:             ", total_sig, "\n")
cat("=========================\n\n")

# ============================================================
# 4. Combine all cell types and sort by p-value
# ============================================================
sorted_gene <- rbindlist(all_results, fill = TRUE)
sorted_gene <- sorted_gene[order(p_value)]
cat("\nTotal gene-celltype pairs (p < 0.05):", nrow(sorted_gene), "\n")

# ============================================================
# 5. Add gene annotation columns
# ============================================================
# Match by gene symbol (CT weights are symbol-named; original matched ENSG on genes$SYM,
# here we match symbol on genes$ID = V4)
transcripts <- sub("\\..*", "", sorted_gene$gene_id)
m <- match(transcripts, genes$ID)

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
  cat("WARNING:", n_na, "gene-celltype pairs could not be matched to gene_annotation.txt.gz\n")
  cat("  These will have NA in CHR/P0/P1/START/STOP columns\n")
}

# ============================================================
# 6. Write output
# ============================================================
write.table(sorted_gene, file = output_file,
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
cat("Written to:", output_file, "\n")
