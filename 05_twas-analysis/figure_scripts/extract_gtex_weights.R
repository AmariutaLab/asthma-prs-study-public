#!/usr/bin/env Rscript
# Extract per-tissue GTEx FUSION eQTL weights for one gene into a long CSV,
# for the eQTL-weight "cards" figure (see build_eqtl_weight_cards.py).
#
# Reads FUSION `.wgt.RDat` models (objects: wgt.matrix [SNP x method], snps,
# cv.performance) from one weight directory per tissue, pulls a single method
# column (default BLUP -- dense, so every SNP has a value), and writes:
#     tissue, snp, weight
#
# Paths/params via environment variables (defaults = local laptop layout):
#     HARTWELL_ROOT   default /Users/nancyh/Desktop/hartwell
#     WEIGHTS_DIR     default $HARTWELL_ROOT/TCSC_asthma/TCSC_weight_files/weights/v8_320EUR
#     ENSG            default ENSG00000073605   (GSDMB, 17q21 asthma locus)
#     METHOD          default blup
#     TISSUES         comma-separated META_ tissue suffixes
#     OUT_CSV         default $TMPDIR/gsdmb_blup4.csv

root    <- Sys.getenv("HARTWELL_ROOT", "/Users/nancyh/Desktop/hartwell")
wdir    <- Sys.getenv("WEIGHTS_DIR",
                      file.path(root, "TCSC_asthma/TCSC_weight_files/weights/v8_320EUR"))
ens     <- Sys.getenv("ENSG", "ENSG00000073605")
meth    <- Sys.getenv("METHOD", "blup")
tissues <- strsplit(Sys.getenv("TISSUES",
                    "Lung,Whole_Blood,Esophagus_Mucosa,Thyroid,Artery_Aorta"), ",")[[1]]
out_csv <- Sys.getenv("OUT_CSV",
                      file.path(Sys.getenv("TMPDIR", "/tmp"), "gsdmb_blup4.csv"))

res <- list()
for (t in tissues) {
  d <- file.path(wdir, paste0("META_", t))
  f <- list.files(d, pattern = paste0("\\.", ens, "\\..*\\.wgt\\.RDat$"), full.names = TRUE)
  if (length(f) == 0) { message("no model for tissue: ", t); next }
  e <- new.env(); load(f[1], envir = e)
  wm <- as.matrix(e$wgt.matrix)
  if (!(meth %in% colnames(wm))) { message("method '", meth, "' absent for ", t); next }
  res[[t]] <- data.frame(tissue = t, snp = rownames(wm),
                         weight = as.numeric(wm[, meth]), stringsAsFactors = FALSE)
}
out <- do.call(rbind, res)
write.csv(out, out_csv, row.names = FALSE)
cat("wrote", nrow(out), "rows across", length(res), "tissues ->", out_csv, "\n")
