library(data.table)

# directory containing CT gene-list files
input_dir <- "${PROJECT_DIR}/TWAS/gene_pt_ct_new/twas_gene_lists"

# list all CT files (e.g., CD14_Mono_gene_list_p5e-05.txt)
# Exclude already-split per-chromosome files (prefix chr*) so re-runs are idempotent
file_list <- list.files(input_dir,
                        pattern = "^[^c].*_gene_list_.*\\.txt$",
                        full.names = TRUE)
file_list <- file_list[!grepl("^chr", basename(file_list))]

for (file in file_list) {

  fname <- basename(file)
  base  <- sub("\\.txt$", "", fname)

  # CT row layout: SYM CHR START STOP ENSG  (col 1 = gene symbol = downstream ID)
  dt <- fread(file, header = FALSE)
  colnames(dt) <- c("SYM", "CHR", "START", "STOP", "ENSG")

  # Split by chromosome
  chrs <- unique(dt$CHR)

  for (chr in chrs) {
    dt_chr <- dt[CHR == chr, ]
    gene_symbols <- dt_chr$SYM

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
