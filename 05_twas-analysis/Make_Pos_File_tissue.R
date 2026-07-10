library(data.table)

# ============================================================
# Build .pos files for GTEx tissue weights (TWAS)
#   Tissue list:    TissueGroups.txt (with N_EUR for the sample-size column)
#   Gene filter:    TranscriptsIn<T>Model_keep.txt (heritable-genes list from
#                   TCSC) intersected with non-NA TWAS.Z in
#                   Marginal_alphas_<trait>_<T>.txt.gz (row-aligned with
#                   TranscriptsIn<T>Model.txt). Both filters together match
#                   the prior local pipeline (the marginal-alphas filter is
#                   what the previous "p_value < 1" expression encoded — it
#                   drops genes whose TWAS regression returned NA).
#   Weights:        /expanse/.../gtex/weights/{v8_320EUR/META_<T>|v8_allEUR_<T>_blup}/<T>.<gene>.wgt.RDat
#   Annotation:     /expanse/.../FOCUS/gene_annotation.txt.gz
#   N per tissue: 320 for every tissue (uniform).
# ============================================================

trait         <- "UKB_460K.disease_ASTHMA_DIAGNOSED"
tissue_groups <- "${REF_DIR}/TCSC/analysis/TissueGroups.txt"
herit_dir     <- "${REF_DIR}/TCSC/weights/heritablegenes"
twas_dir      <- "${REF_DIR}/TCSC/twas_statistics"
annot_file    <- "${PROJECT_DIR}/FOCUS/gene_annotation.txt.gz"
outdir        <- "${PROJECT_DIR}/TWAS/TWAS_pos"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Tissue list + per-tissue sample size
# ------------------------------------------------------------
y          <- fread(tissue_groups, header = TRUE)
tissues    <- unique(y$MetaTissue)
n_eqtl     <- sapply(tissues, function(t) sum(y$N_EUR[y$MetaTissue == t]))
small_flag <- n_eqtl < 320
cat("Found", length(tissues), "tissues (",
    sum(small_flag), "small -> Nall blup model)\n")

# ------------------------------------------------------------
# Gene annotation: V1=chr, V2=start, V3=stop, V4=symbol, V7=ENSG (versioned)
# ------------------------------------------------------------
genes_full <- fread(annot_file, sep = "\t", header = FALSE)
genes <- genes_full[, .(V1, V2, V3, V4, V7)]
setnames(genes, c("CHR", "START", "STOP", "SYM", "ENSG"))
genes[, ENSG := sub("\\..*", "", ENSG)]
genes <- unique(genes, by = "ENSG")

# ------------------------------------------------------------
# Iterate tissues
# ------------------------------------------------------------
for (i in seq_along(tissues)) {
  tissue <- tissues[i]

  if (small_flag[i]) {
    prefix         <- "Nall"
    twas_subfolder <- "allEUR_GTEx"
    weights_subdir <- paste0("v8_allEUR_", tissue, "_blup/", tissue, ".")
  } else {
    prefix         <- "N320"
    twas_subfolder <- "320EUR_GTEx"
    weights_subdir <- paste0("v8_320EUR/META_", tissue, "/", tissue, ".")
  }
  n_t <- 320L

  transcripts_file <- file.path(herit_dir, prefix,
                                paste0("TranscriptsIn", tissue, "Model.txt"))
  keep_file        <- file.path(herit_dir, prefix,
                                paste0("TranscriptsIn", tissue, "Model_keep.txt"))
  twas_file        <- file.path(twas_dir, twas_subfolder, trait,
                                paste0("Marginal_alphas_", trait, "_", tissue, ".txt.gz"))

  if (!file.exists(keep_file) || !file.exists(transcripts_file) || !file.exists(twas_file)) {
    message("WARNING: missing input for ", tissue, " — skipping"); next
  }

  # Gene-level z-scores: full transcripts list row-aligned with marginal alphas
  full   <- fread(transcripts_file, header = FALSE)   # column V1 = gene_id (ENSG.version)
  twas_z <- fread(twas_file,        header = FALSE)   # column V1 = TWAS z-score
  setnames(full, "V1", "gene_id")
  full[, twas_z := twas_z$V1]

  # Intersect with heritable-genes keep list + drop NA z (matches prior
  # local code's `p_value < 1` filter, which excluded NA marginal alphas)
  keep <- fread(keep_file, header = FALSE)
  full <- full[gene_id %in% keep$V1 & !is.na(twas_z)]

  gene_ids      <- full$gene_id
  ensg_stripped <- sub("\\..*", "", gene_ids)
  m             <- match(ensg_stripped, genes$ENSG)

  dt <- data.table(
    PANEL = tolower(tissue),
    WGT   = paste0(weights_subdir, gene_ids, ".wgt.RDat"),
    ID    = genes$SYM[m],
    CHR   = genes$CHR[m],
    P0    = genes$START[m] - 1000000,    # 1Mb buffer (annotation may differ from plink build)
    P1    = genes$STOP[m]  + 1000000,
    N     = n_t
  )

  # Strip any trailing ".N" suffix on the gene symbol (matches prior local code)
  dt[, ID  := sub("\\..*", "", ID)]
  dt[, CHR := suppressWarnings(as.numeric(gsub("chr", "", CHR)))]

  missing <- dt[is.na(CHR)]
  if (nrow(missing) > 0) {
    cat(sprintf("  %s: %d / %d genes unmatched in annotation (dropped)\n",
                tissue, nrow(missing), nrow(dt)))
  }
  dt <- dt[!is.na(CHR) & CHR > 0]

  setorder(dt, CHR, P0)
  outfile <- file.path(outdir, paste0("twas_", tolower(tissue), ".pos"))
  fwrite(dt, file = outfile, sep = "\t", col.names = TRUE)

  cat(sprintf("%-40s  %5d genes -> %s (N=%d)\n",
              tissue, nrow(dt), outfile, n_t))
}

cat("\nDone. .pos files in:", outdir, "\n")
