# ============================================================
# Rebuild the 39-tissue TWAS path table consumed by stage 1 of the PTRS
# construction step (08_ptrs-construction/ptrs_score_v4.sbatch). For every
# tissue in TissueGroups.txt this script:
#   1. Loads the per-tissue TWAS .dat files from TWAS_results/
#   2. Writes a row-aligned pair of files (one per tissue) into rebuilt/:
#        - rebuilt/gtex_twas/Marginal_alphas_NEW_TWAS_<Tissue>.txt.gz  (TWAS.Z column)
#        - rebuilt/transcripts/TranscriptsIn<Tissue>Model.txt          (gene_id column)
#   3. Emits the 6-column path table at rebuilt/twas_path_table_39_new.csv
#      with columns: tissue, gtex_twas, transcripts, transcript_keep_file,
#      weights, PANEL. transcript_keep_file points at the heritable-genes
#      filter list from the upstream TCSC pipeline (N320/Nall by tissue size).
# ============================================================

library(data.table)

# ============================================================
# Paths
# ============================================================
twas_dir       <- "${PROJECT_DIR}/TWAS/TWAS_results"
tissue_file    <- "${REF_DIR}/TCSC/analysis/TissueGroups.txt"
herit_dir      <- "${REF_DIR}/TCSC/weights/heritablegenes"

# Output directories for the rebuilt files
out_base       <- "${PROJECT_DIR}/TWAS/rebuilt"
out_twas_dir   <- file.path(out_base, "gtex_twas")     # holds Marginal_alphas_*.txt.gz
out_trans_dir  <- file.path(out_base, "transcripts")   # holds TranscriptsIn*Model.txt
out_table      <- file.path(out_base, "twas_path_table_39_new.csv")

dir.create(out_twas_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_trans_dir, recursive = TRUE, showWarnings = FALSE)

# Original GTEx weights base path (kept the same — your weights live here)
weights_base_320  <- "${REF_DIR}/gtex/weights/v8_320EUR/META_"
weights_base_nall <- "${REF_DIR}/gtex/weights/v8_allEUR_"

# ============================================================
# 1. Tissue groups
# ============================================================
y <- fread(tissue_file, header = TRUE)
tissues    <- unique(y$MetaTissue)
n_eqtl     <- sapply(tissues, function(t) sum(y$N_EUR[y$MetaTissue == t]))
small_flag <- n_eqtl < 320

cat("Found", length(tissues), "unique MetaTissues\n")

# ============================================================
# 2. For each tissue: load TWAS, write gtex_twas + transcripts
# ============================================================
table_rows <- list()

for (i in seq_along(tissues)) {
  tissue       <- tissues[i]
  tissue_lower <- tolower(tissue)

  # Load TWAS data (prefer _all.dat, else combine _chr*.dat)
  all_file <- file.path(twas_dir, paste0("twas_", tissue_lower, "_all.dat"))
  if (file.exists(all_file)) {
    dat <- fread(all_file, header = TRUE)
  } else {
    chr_files <- list.files(twas_dir,
                            pattern = paste0("^twas_", tissue_lower, "_chr[0-9]+\\.dat$"),
                            full.names = TRUE)
    if (length(chr_files) == 0) {
      message("WARNING: no TWAS files for tissue: ", tissue)
      next
    }
    dat <- rbindlist(lapply(chr_files, fread, header = TRUE))
  }

  # Extract ENSG ID (with version) from the FILE column
  dat[, gene_id := sub("^.*\\.(ENSG[0-9]+\\.[0-9]+)\\.wgt\\.RDat$", "\\1", FILE)]

  # Drop any rows with missing TWAS.Z
  dat <- dat[!is.na(TWAS.Z)]

  # ----------------------------------------------------------
  # Write parallel files (matching original format):
  #   gtex_twas: gzipped, one z-score per line, no header
  #   transcripts: one ENSG ID per line, no header
  # The two files are row-aligned (cbind-able).
  # ----------------------------------------------------------
  prefix <- ifelse(small_flag[i], "Nall", "N320")

  gtex_twas_file   <- file.path(out_twas_dir,
                                paste0("Marginal_alphas_NEW_TWAS_", tissue, ".txt.gz"))
  transcripts_file <- file.path(out_trans_dir,
                                paste0("TranscriptsIn", tissue, "Model.txt"))

  fwrite(dat[, .(TWAS.Z)],
         file = gtex_twas_file,
         col.names = FALSE, quote = FALSE, sep = "\t")
  fwrite(dat[, .(gene_id)],
         file = transcripts_file,
         col.names = FALSE, quote = FALSE, sep = "\t")

  # Use the original keep file (heritability filter list — unchanged)
  transcript_keep_file <- file.path(herit_dir, prefix,
                                    paste0("TranscriptsIn", tissue, "Model_keep.txt"))

  # Weights base path
  if (small_flag[i]) {
    weights_path <- paste0(weights_base_nall, tissue, "_blup/", tissue, ".")
  } else {
    weights_path <- paste0(weights_base_320, tissue, "/", tissue, ".")
  }

  table_rows[[length(table_rows) + 1]] <- data.table(
    tissue               = tissue,
    gtex_twas            = gtex_twas_file,
    transcripts          = transcripts_file,
    transcript_keep_file = transcript_keep_file,
    weights              = weights_path,
    PANEL                = tissue_lower
  )

  cat(sprintf("%-40s  %d genes written\n", tissue, nrow(dat)))
}

# ============================================================
# 3. Write twas_path_table_39.csv
# ============================================================
table_dt <- rbindlist(table_rows)
fwrite(table_dt, file = out_table, quote = TRUE)
cat("\nWrote path table:", out_table, "\n")
cat("Total tissues:", nrow(table_dt), "\n")

