#!/bin/bash
# =============================================================================
# run_demo.sh — self-contained smoke test for TWAS p-value LD-clumping (tissue)
# =============================================================================
#
# Reproduces the core method of ld_clumping_tissue.R + split_tissue.R on the
# committed demo fixture (Whole_Blood, 12 chr1 genes), so anyone can confirm the
# gene-selection logic runs without the full GTEx predicted-expression release.
#
# Method (identical to ld_clumping_tissue.R):
#   1. gene-gene correlation from the predicted-expression matrix (cor(df1))
#   2. two-sided p from TWAS Z  (p = 2*pnorm(|Z|, lower.tail=FALSE))
#   3. sort genes by p ascending
#   4. greedy LD clump: keep a gene iff r^2 < 0.1 vs every already-kept gene
#   5. intersect with the heritability keep-list
#   6. emit gene lists at p in {5e-5, 5e-4, 5e-3, 0.05}, then split by chromosome
#
# Tissue track only. The 17-cell-type (OneK1K) track is PENDING (see README).
#
# DEMO DATA (see sample_data/README.md):
#   designmat_Whole_Blood_demo.RData             REAL GTEx v8 predicted expression
#       (489 samples x 12 genes) — object `df1`.
#   Marginal_alphas_NEW_TWAS_Whole_Blood.demo.txt.gz   TWAS Z from 05's demo run.
#   TranscriptsInWhole_BloodModel.demo.txt       gene ids row-aligned with df1.
#   TranscriptsInWhole_BloodModel_keep.demo.txt  heritability keep-list.
#   gene annotation is reused from ../07_focus-finemapping/data/gene_annotation.txt.gz
#
# REQUIREMENTS: R with `data.table`.
# =============================================================================

set -e -o pipefail
RSCRIPT="${RSCRIPT:-Rscript}"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "=========================================================================="
echo " TWAS P+T LD-clumping — demo smoke test (Whole_Blood, 12 chr1 genes)"
echo "   Rscript : $(command -v "${RSCRIPT}" || echo NOT FOUND)"
echo "=========================================================================="

"${RSCRIPT}" - "${SCRIPT_DIR}" <<'RS'
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE); DIR <- args[1]
SD  <- file.path(DIR, "sample_data")
OUT <- file.path(DIR, "demo_run", "twas_gene_lists")
unlink(file.path(DIR, "demo_run"), recursive = TRUE); dir.create(OUT, recursive = TRUE)

anno_file <- file.path(DIR, "..", "07_focus-finemapping", "data", "gene_annotation.txt.gz")
genes <- fread(anno_file, sep = "\t", header = FALSE)[, .(CHR=V1, START=V2, STOP=V3, ID=V4, SYM=V7)]
genes[, SYM := sub("\\..*$", "", SYM)][, CHR := as.character(CHR)]

# --- load fixture ---
load(file.path(SD, "designmat_Whole_Blood_demo.RData"))            # df1
tr   <- fread(file.path(SD, "TranscriptsInWhole_BloodModel.demo.txt"), header = FALSE)$V1
keep <- sub("\\..*$", "", fread(file.path(SD, "TranscriptsInWhole_BloodModel_keep.demo.txt"), header = FALSE)$V1)
z    <- fread(file.path(SD, "Marginal_alphas_NEW_TWAS_Whole_Blood.demo.txt.gz"), header = FALSE)$V1
base <- sub("\\..*$", "", tr)

# --- 1. correlation matrix, labelled by gene ---
corr <- cor(df1); rownames(corr) <- base; colnames(corr) <- base

# --- 2-3. p from Z, sort ascending ---
d <- data.table(gene_id = base, Z = z)[, p_value := 2 * pnorm(abs(Z), lower.tail = FALSE)]
d <- d[order(p_value)]

# --- 4. greedy LD clump (r^2 < 0.1) ---
sel <- character(0)
for (g in d$gene_id) {
  if (length(sel) == 0) { sel <- g; next }
  if (all((corr[g, sel])^2 < 0.1, na.rm = TRUE)) sel <- c(sel, g)
}
# --- 5. heritability keep-list ---
sel <- sel[sel %in% keep]
cat("Clumped genes kept:", length(sel), "of", nrow(d), "\n\n")

# --- 6. threshold + split by chr ---
for (thr in c(5e-5, 5e-4, 5e-3, 0.05)) {
  gl <- d[gene_id %in% sel & p_value < thr]$gene_id
  m  <- merge(data.table(gene_id = gl), genes, by.x = "gene_id", by.y = "SYM", all.x = TRUE)
  lab <- gsub("\\.", "_", as.character(thr))
  fwrite(m, file.path(OUT, paste0("Whole_Blood_gene_list_p", lab, ".txt")),
         sep = " ", quote = FALSE, col.names = FALSE)
  for (chr in unique(m$CHR[!is.na(m$CHR)]))
    fwrite(m[CHR == chr, .(ID)], file.path(OUT, paste0(chr, "_Whole_Blood_gene_list_p", lab, ".txt")),
           col.names = FALSE, quote = FALSE)
  cat(sprintf("p<%-6s : %d genes -> Whole_Blood_gene_list_p%s.txt\n", thr, length(gl), lab))
}
cat("\nOutputs in", OUT, "\n")
RS

echo; echo "Demo complete."
