#!/usr/bin/env Rscript
# Extract per-cell-type OneK1K FUSION eQTL weights for one gene into a long CSV,
# for the single-cell eQTL-weight "cards" figure (build_eqtl_weight_cards_1k1k.py).
#
# Reads FUSION `.wgt.RDat` models (objects: wgt.matrix [SNP x method], snps,
# cv.performance) from one weight directory per OneK1K cell type, laid out as
#     <MAFOCUS_DIR>/<CellType>/<GENE>.wgt.RDat     (weights keyed by gene SYMBOL)
# pulls a single method column (default BLUP -- dense, so every SNP has a value),
# and writes:  celltype, snp, weight
#
# Paths/params via environment variables (defaults = local laptop layout):
#     HARTWELL_ROOT   default /Users/nancyh/Desktop/hartwell
#     MAFOCUS_DIR     default $HARTWELL_ROOT/FOCUS/ma-focus
#     GENE            default HLA-DPB1   (cell-type panel focal gene; SYMBOL)
#     METHOD          default blup
#     CELLTYPES       comma-separated CamelCase cell-type dir names
#     OUT_CSV         default $TMPDIR/hla_dpb1_onek1k_blup.csv

root <- Sys.getenv("HARTWELL_ROOT", "/Users/nancyh/Desktop/hartwell")
wdir <- Sys.getenv("MAFOCUS_DIR", file.path(root, "FOCUS/ma-focus"))
gene <- Sys.getenv("GENE", "HLA-DPB1")
meth <- Sys.getenv("METHOD", "blup")
cts  <- strsplit(Sys.getenv("CELLTYPES",
                 "B_memory,CD4_CTL,CD4_TCM,NK"), ",")[[1]]
out_csv <- Sys.getenv("OUT_CSV",
                 file.path(Sys.getenv("TMPDIR", "/tmp"), "hla_dpb1_onek1k_blup.csv"))

res <- list()
for (ct in cts) {
  f <- file.path(wdir, ct, paste0(gene, ".wgt.RDat"))
  if (!file.exists(f)) { message("no model for cell type: ", ct); next }
  e <- new.env(); load(f, envir = e)
  wm <- as.matrix(e$wgt.matrix)
  if (!(meth %in% colnames(wm))) { message("method '", meth, "' absent for ", ct); next }
  res[[ct]] <- data.frame(celltype = ct, snp = rownames(wm),
                          weight = as.numeric(wm[, meth]), stringsAsFactors = FALSE)
}
out <- do.call(rbind, res)
write.csv(out, out_csv, row.names = FALSE)
cat("wrote", nrow(out), "rows across", length(res), "cell types ->", out_csv, "\n")
