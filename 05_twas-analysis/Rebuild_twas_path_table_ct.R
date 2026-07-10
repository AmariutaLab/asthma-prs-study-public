# ============================================================
# Rebuild the 17-cell-type TWAS path table consumed by stage 1 of the PTRS
# construction step (08_ptrs-construction/ptrs_score_v5.sbatch). For every
# cell type with a .pos file in TWAS_pos_ct/ this script:
#   1. Loads the per-CT TWAS .dat files from TWAS_results_ct/
#   2. Filters down to the gene-CT pairs that passed QC_twas_ct.R
#      (HSQ>0, nonzero predictor variance, MODELCV.R2>0)
#   3. Writes a row-aligned pair of files (one per CT) into rebuilt_ct/:
#        - rebuilt_ct/ct_twas/Marginal_alphas_NEW_TWAS_<ct>.txt.gz  (TWAS.Z column)
#        - rebuilt_ct/transcripts/TranscriptsIn<ct>Model.txt        (gene_id column)
#   4. Emits the 5-column path table at rebuilt_ct/twas_path_table_ct_new.csv
#      with columns: tissue, gtex_twas, transcripts, weights, PANEL.
#      (No transcript_keep_file column: QC has already been applied.)
# ============================================================

library(data.table)

# ============================================================
# Paths
# ============================================================
twas_dir       <- "${PROJECT_DIR}/TWAS/TWAS_results_ct"
pos_dir        <- "${PROJECT_DIR}/TWAS/TWAS_pos_ct"
qc_genelist    <- "${PROJECT_DIR}/TWAS/twas_qc_genelist_ct.csv"

# Output directories for the rebuilt files
out_base       <- "${PROJECT_DIR}/TWAS/rebuilt_ct"
out_twas_dir   <- file.path(out_base, "ct_twas")       # holds Marginal_alphas_*.txt.gz
out_trans_dir  <- file.path(out_base, "transcripts")   # holds TranscriptsIn*Model.txt
out_table      <- file.path(out_base, "twas_path_table_ct_new.csv")

dir.create(out_twas_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_trans_dir, recursive = TRUE, showWarnings = FALSE)

# Cell-type weights base path (FUSION weights live here)
weights_base <- "${REF_DIR}/1k1k_cell_type/celltype_weights_v2/"

# ============================================================
# 1. Cell types (derived from pos files)
# ============================================================
pos_files <- list.files(pos_dir, pattern = "^twas_.*\\.pos$", full.names = TRUE)
tissues   <- gsub("^twas_|\\.pos$", "", basename(pos_files))

cat("Found", length(tissues), "cell types\n")

# Load QC gene list (PANEL, gene_id) produced by QC_twas_ct.R
qc <- fread(qc_genelist)
qc[, key := paste0(PANEL, "|", gene_id)]
cat("QC gene list loaded:", nrow(qc), "gene-CT pairs\n")

# ============================================================
# 2. For each cell type: load TWAS, write gtex_twas + transcripts
# ============================================================
table_rows <- list()

for (i in seq_along(tissues)) {
  tissue       <- tissues[i]
  tissue_lower <- tissue

  # Recover case-preserved CT directory name from pos file's first WGT entry
  pos    <- fread(file.path(pos_dir, paste0("twas_", tissue_lower, ".pos")), header = TRUE)
  ct_dir <- sub("/.*$", "", pos$WGT[1])

  # Load TWAS data (prefer _all.dat, else combine _chr*.dat)
  all_file <- file.path(twas_dir, paste0("TWAS_metaGWAS_", tissue_lower, "_all.dat"))
  if (file.exists(all_file)) {
    dat <- fread(all_file, header = TRUE)
  } else {
    chr_files <- list.files(twas_dir,
                            pattern = paste0("^TWAS_metaGWAS_", tissue_lower, "_chr[0-9]+\\.dat$"),
                            full.names = TRUE)
    if (length(chr_files) == 0) {
      message("WARNING: no TWAS files for cell type: ", tissue)
      next
    }
    dat <- rbindlist(lapply(chr_files, fread, header = TRUE))
  }

  # Gene ID = symbol (CT weights are symbol-named; FUSION propagates ID column from pos)
  dat[, gene_id := ID]

  # Keep only gene-CT pairs that passed QC (HSQ>0, nonzero var, MODELCV.R2>0)
  dat <- dat[paste0(tissue_lower, "|", gene_id) %in% qc$key]

  # ----------------------------------------------------------
  # Write parallel files (matching original format):
  #   gtex_twas: gzipped, one z-score per line, no header
  #   transcripts: one gene ID per line, no header
  # The two files are row-aligned (cbind-able).
  # ----------------------------------------------------------
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

  # Weights base path
  weights_path <- paste0(weights_base, ct_dir, "/")

  table_rows[[length(table_rows) + 1]] <- data.table(
    tissue               = tissue,
    gtex_twas            = gtex_twas_file,
    transcripts          = transcripts_file,
    weights              = weights_path,
    PANEL                = tissue_lower
  )

  cat(sprintf("%-40s  %d genes written\n", tissue, nrow(dat)))
}

# ============================================================
# 3. Write twas_path_table_ct.csv
# ============================================================
table_dt <- rbindlist(table_rows)
fwrite(table_dt, file = out_table, quote = TRUE)
cat("\nWrote path table:", out_table, "\n")
cat("Total cell types:", nrow(table_dt), "\n")
