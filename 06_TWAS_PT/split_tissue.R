library(data.table)

# directory containing tissue gene-list files
input_dir <- "${PROJECT_DIR}/TWAS/gene_pt_new/twas_gene_lists"

# list all tissue files (e.g., Brain_Cerebellum_gene_list_p5e-05.txt)
# Exclude already-split per-chromosome files (prefix chr*) so re-runs are idempotent
file_list <- list.files(input_dir,
                        pattern = "^[^c].*_gene_list_.*\\.txt$",
                        full.names = TRUE)
file_list <- file_list[!grepl("^chr", basename(file_list))]

for (file in file_list) {

  # Extract the base filename (e.g. Brain_Cerebellum_gene_list_p5e-05.txt)
  fname <- basename(file)

  # Extract tissue + threshold part from filename
  # e.g. "Brain_Cerebellum_gene_list_p5e-05"
  base <- sub("\\.txt$", "", fname)

  # Load data
  dt <- fread(file, header = FALSE)
  colnames(dt) <- c("SYM", "CHR", "START", "STOP", "gene_id")

  # Split by chromosome
  chrs <- unique(dt$CHR)

  for (chr in chrs) {
    dt_chr <- dt[CHR == chr, ]

    # Keep only SYM
    gene_symbols <- dt_chr$SYM

    # Output file name:
    # chr16_Brain_Cerebellum_gene_list_p5e-05.txt
    out_file <- file.path(
      input_dir,
      paste0(chr, "_", base, ".txt")
    )

    fwrite(
      data.table(SYM = gene_symbols),
      out_file,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    cat("Created:", out_file, " (", length(gene_symbols), "genes )\n")
  }
}
