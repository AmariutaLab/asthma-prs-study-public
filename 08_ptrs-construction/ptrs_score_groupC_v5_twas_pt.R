#!/usr/bin/env Rscript
################################################################################
# ptrs_score_groupC_v5_twas_pt.R
#
# Per-chromosome PTRS scoring for one cell type (17 OneK1K CT variant),
# TWAS P+T (p-value-thresholding) variant, called by
# `ptrs_score_v5_twas_pt.sbatch`. Same scoring loop as `ptrs_score_groupC_v5.R`
# (the FOCUS-genelist variant); the difference is the gene-list source -- here
# the `gene_list_template` argument points at the LD-clumped, per-pval TWAS P+T
# gene lists under TWAS/gene_pt_ct_new/twas_gene_lists/, whose filenames use
# the case-preserved cell-type name. For each significant gene in the
# chromosome, this script:
#   1. Loads the FUSION/groupC-trained weight file for the gene
#   2. Picks the best model (top1 / lasso / enet / blup) by cv.performance
#      p-value -- except for groupC tissues whose weights have a single column,
#      where that column's beta is used directly
#   3. Strand-flips / allele-swaps SNPs in the cohort BIM to match the
#      effect-allele in the weight file, writes the updated BIM in place
#   4. Runs `plink --score` to get one per-sample SCORE per gene
#   5. Sums those scores across genes, with three weighting variants:
#        - all_data        : raw sum
#        - all_data_keep   : sum weighted by signed gene z-score
#        - all_data_sign   : sum weighted by sign(gene z-score) only
#   6. Writes <output_dir>/chr<chr>_<tissue>_all_data{,_keep,_sign}.rds
#
# Author: Nan Huang (n5huang@ucsd.edu)
#
# Usage:
#   Rscript ptrs_score_groupC_v5_twas_pt.R \
#       <chr> <transcript_file> <gtex_file> <transcript_filter_file> \
#       <weight_file_basename> <output_dir> <tissue> <pval> \
#       <plink_basename> <gene_list_template>
#
# Arguments:
#   chr                     Chromosome number (1-22)
#   transcript_file         Per-tissue transcript list (one gene_id per line)
#   gtex_file               Per-tissue TWAS z-score column (same order as transcript_file)
#   transcript_filter_file  Subset of gene_ids to retain (pass = transcript_file for no-op)
#   weight_file_basename    Prefix for per-gene .wgt.RDat files; script appends "<gene_id>.wgt.RDat"
#   output_dir              Directory for per-gene .profile + per-chr .rds outputs
#   tissue                  Tissue / cell-type label (used in gene-list path + output filenames)
#   pval                    P-value tag (used in gene-list template substitution)
#   plink_basename          PLINK bfile prefix; script appends "<chr>" to get bim/bed/fam
#   gene_list_template      Path to per-(chr, cell-type, pval) gene list with
#                           placeholders {chr}, {tissue}, {pval}. {tissue} is
#                           filled with the case-preserved cell-type name
#                           (weights-dir basename, e.g. "B_intermediate"), not
#                           the lowercase arg-7 label. Example:
#                             .../TWAS/gene_pt_ct_new/twas_gene_lists/chr{chr}_{tissue}_gene_list_p{pval}.txt
################################################################################

args <- commandArgs(trailingOnly = TRUE)

EXPECTED_ARGS <- 10
write(paste("[USER INPUT] Number of Arguments:", length(args)), stdout())
if (length(args) != EXPECTED_ARGS) {
    if (length(args) > 0) {
        write("[USER INPUT] Arguments provided:", stdout())
        for (i in seq_along(args)) {
            write(paste0("  Arg", i, ": ", args[[i]]), stdout())
        }
    }
    stop(sprintf("[ERROR] Expected %d arguments, got %d", EXPECTED_ARGS, length(args)),
         call. = FALSE)
}

input_chrm        <- args[1]
input_tran        <- args[2]
input_gtex        <- args[3]
input_filt        <- args[4]
input_wght        <- args[5]
output_dir        <- args[6]
tissue_type       <- args[7]
pval              <- args[8]
plink_basename    <- args[9]
gene_list_template <- args[10]

write(paste("[USER INPUT] Chromosome:", input_chrm), stdout())
write(paste("[USER INPUT] Transcript File:", input_tran), stdout())
write(paste("[USER INPUT] GTEX File:", input_gtex), stdout())
write(paste("[USER INPUT] Transcript Filter File:", input_filt), stdout())
write(paste("[USER INPUT] Weight File Basename:", input_wght), stdout())
write(paste("[USER INPUT] Output Location:", output_dir), stdout())
write(paste("[USER INPUT] Tissue:", tissue_type), stdout())
write(paste("[USER INPUT] P value:", pval), stdout())
write(paste("[USER INPUT] PLINK basename:", plink_basename), stdout())
write(paste("[USER INPUT] Gene-list template:", gene_list_template), stdout())

################################################################################
suppressPackageStartupMessages({
    library(data.table)
    library(glmnet)
    library(dplyr)
})

# Tissues whose weight matrix has a single column (use that column directly)
GROUP_C_TISSUES <- c(
    "Brain_Limbic", "Brain_Cortex", "Brain_BasalGanglia", "Brain_Cerebellum",
    "Colon_Sigmoid", "Colon_Transverse"
)

# Substitute {chr}, {tissue}, {pval} in the gene-list template
resolve_template <- function(tpl, chr, tissue, pval) {
    out <- gsub("\\{chr\\}",    chr,    tpl)
    out <- gsub("\\{tissue\\}", tissue, out)
    out <- gsub("\\{pval\\}",   pval,   out)
    out
}

flip_allele <- function(allele) {
    allele <- toupper(allele)
    switch(allele, "A" = "T", "T" = "A", "C" = "G", "G" = "C", allele)
}

################################################################################
# Load + filter transcripts -- once per chromosome
transcripts <- fread(input_tran, header = FALSE)
setnames(transcripts, "V1", "gene_id")
gtex_twas <- fread(input_gtex)
gtex_raw  <- cbind(transcripts, gtex_twas)
write(paste0("[ ", Sys.time(), " ] Loaded Transcript and GTEX Files into [ gtex_raw ]"), stdout())

transcripts_keep <- fread(input_filt, header = FALSE)
gtex_combine <- merge(gtex_raw, transcripts_keep, by.x = "gene_id", by.y = "V1", all.y = TRUE)
gtex_combine <- gtex_combine[complete.cases(gtex_combine), ]
write(paste0("[ ", Sys.time(), " ] Filtered into [ gtex_combine ]"), stdout())

# Load the per-(chr, cell-type) significant-gene list.
# CT P+T gene-list files are keyed by the case-preserved cell-type name (the
# weights-dir basename, e.g. "B_intermediate"); arg-7 `tissue_type` is the
# lowercase label used only for output filenames -- so resolve {tissue} with
# the cased name here.
ct_name_cased  <- basename(sub("/$", "", input_wght))
gene_list_path <- resolve_template(gene_list_template, input_chrm, ct_name_cased, pval)
if (!file.exists(gene_list_path)) {
    stop(paste("Gene list file not found:", gene_list_path), call. = FALSE)
}
gene_lst <- fread(gene_list_path, header = FALSE)
setnames(gene_lst, "V1")
write(paste("Loaded significant genes from:", gene_list_path), stdout())

# Match on version-stripped gene IDs (gene lists are unversioned)
gtex_combine <- gtex_combine %>%
    mutate(gene_id_strip = sub("\\..*", "", gene_id)) %>%
    filter(gene_id_strip %in% gene_lst$V1)

num_significant_genes <- nrow(gtex_combine)
write(paste("Number of significant genes retained:", num_significant_genes), stdout())
if (num_significant_genes == 0) {
    stop("No significant genes found after filtering. Please check the gene list.",
         call. = FALSE)
}

################################################################################
# Resolve cohort genotype files from PLINK basename + chromosome
plink_file    <- paste0(plink_basename, input_chrm)
bim_file_path <- paste0(plink_file, ".bim")
write(paste("[RESOLVED] .bim File Path:", bim_file_path), stdout())
write(paste("[RESOLVED] PLINK File Path:", plink_file), stdout())

################################################################################
all_data       <- NULL
all_data_keep  <- NULL
all_data_sign  <- NULL

for (i in 1:nrow(gtex_combine)) {
    gene_id <- gtex_combine$gene_id[i]
    print(gene_id)

    file_path <- paste0(input_wght, gene_id, ".wgt.RDat")
    if (!file.exists(file_path)) {
        write(paste("File not found for gene:", gene_id), stdout())
        next
    }

    load(file_path)

    bim_file <- fread(bim_file_path, header = FALSE)
    matching_snps <- snps$V2 %in% bim_file$V2
    if (sum(matching_snps) == 0) {
        print(1)
        next
    }

    bim_snps  <- bim_file$V2
    snp_names <- intersect(snps$V2, bim_snps)

    total_snps_snps  <- length(snps$V2)
    total_snps_bim   <- length(bim_snps)
    num_intersect    <- length(snp_names)
    pct_snps_in_snps <- 100 * num_intersect / total_snps_snps
    pct_snps_in_bim  <- 100 * num_intersect / total_snps_bim
    write(sprintf("[ %s ] SNP overlap: %d matched | %.2f%% of snps$V2 | %.2f%% of bim_file$V2",
                  Sys.time(), num_intersect, pct_snps_in_snps, pct_snps_in_bim), stdout())

    if (!(tissue_type %in% GROUP_C_TISSUES)) {
        # ---- Standard 4-model selection ----------------------------------------
        cv_row      <- cv.performance["pval", ]
        pvalues     <- c(cv_row["top1"], cv_row["lasso"], cv_row["enet"], cv_row["blup"])
        model_names <- c("top1", "lasso", "enet", "blup")
        pvalues[is.na(pvalues)] <- -Inf
        min_index <- which.min(pvalues)

        if (min(pvalues) >= 0.05) next  # no model passes p < 0.05

        model <- model_names[min_index]
        print(model)

        beta <- as.vector(wgt.matrix[snp_names, model])

        if (model == "top1") {
            beta_squared   <- beta^2
            max_snp_index  <- which.max(beta_squared)
            beta_filtered  <- beta[max_snp_index]
            snp_name       <- snp_names[max_snp_index]
            matching_rows  <- snps %>% filter(V2 %in% snp_name)
            if (nrow(matching_rows) == 0) {
                cat("No matching SNP found for top1 model:", snp_name, "\n")
            }
            Asthma_score_file <- cbind(matching_rows$V2, matching_rows$V5, beta_filtered)
        } else {
            matching_rows     <- snps %>% filter(V2 %in% snp_names)
            Asthma_score_file <- cbind(matching_rows$V2, matching_rows$V5, beta)
            Asthma_score_file <- Asthma_score_file[!duplicated(Asthma_score_file[, 1]), ]
        }
    } else {
        # ---- groupC tissues: single-column wgt.matrix --------------------------
        beta              <- as.vector(wgt.matrix[snp_names, ])
        matching_rows     <- snps %>% filter(V2 %in% snp_names)
        Asthma_score_file <- cbind(matching_rows$V2, matching_rows$V5, beta)
        Asthma_score_file <- Asthma_score_file[!duplicated(Asthma_score_file[, 1]), ]
    }

    if (nrow(Asthma_score_file) <= 0) {
        print("no valid snps for score file")
        next
    }
    if (any(is.na(Asthma_score_file[, 3]))) {
        Asthma_score_file <- Asthma_score_file[!is.na(Asthma_score_file[, 3]), ]
    }
    if (nrow(Asthma_score_file) <= 0) {
        print("no valid snps after NA filter")
        next
    }

    Asthma_score_file_path <- paste0(output_dir, "/chr", input_chrm,
                                     "_Asthma_score_file_", gene_id, ".txt")
    write.table(Asthma_score_file, file = Asthma_score_file_path,
                row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
    print(head(Asthma_score_file))

    # ------------------------------------------------------------------ strand fix
    score <- fread(Asthma_score_file_path, header = FALSE)
    colnames(score) <- c("SNP", "EFFECT_GTEX", "BETA")
    # A1 = ALT (effect), A2 = REF in our cohort BIMs
    colnames(bim_file) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

    merged <- merge(bim_file, score, by = "SNP")
    diff_both_count <- sum(merged$A1 != merged$EFFECT_GTEX &
                           merged$A2 != merged$EFFECT_GTEX, na.rm = TRUE)
    cat("Number of SNPs where both A1 and A2 != EFFECT_GTEX:", diff_both_count, "\n")

    mismatch_idx <- !is.na(merged$EFFECT_GTEX) & merged$A1 != merged$EFFECT_GTEX

    # Case 1: A2 == EFFECT_GTEX -> swap A1/A2
    swap_idx <- mismatch_idx & merged$A2 == merged$EFFECT_GTEX
    if (any(swap_idx)) {
        merged[swap_idx, c("A1", "A2") := .(A2, A1)]
    }
    # Case 2: otherwise strand-flip both
    strand_flip_idx <- mismatch_idx & !swap_idx
    if (any(strand_flip_idx)) {
        merged[strand_flip_idx, A1 := EFFECT_GTEX]
        merged[strand_flip_idx, A2 := sapply(EFFECT_GTEX, flip_allele)]
    }

    bim_updated <- merged[, .(CHR, SNP, CM, BP, A1, A2)]
    merged2     <- merge(bim_updated, score, by = "SNP")
    diff_both_count2 <- sum(merged2$A2 != merged2$EFFECT_GTEX &
                            merged2$A1 != merged2$EFFECT_GTEX, na.rm = TRUE)
    cat("Number of SNPs where both A2 and A1 != EFFECT_GTEX:", diff_both_count2, "\n")

    # Update only overlapping SNPs in the cohort BIM and write it back in place
    setkey(bim_file, SNP)
    setkey(bim_updated, SNP)
    bim_file[bim_updated, `:=`(A1 = i.A1, A2 = i.A2)]
    fwrite(bim_file, bim_file_path, sep = "\t", col.names = FALSE)

    # ----------------------------------------------------------------- PLINK score
    plink_command <- paste(
        "${PROJECT_DIR}/./plink --bfile", plink_file,
        "--out", paste0(output_dir, "/", gene_id, "_out"),
        "--score", Asthma_score_file_path,
        "1 2 3"
    )
    system(plink_command)

    gene_data <- fread(paste0(output_dir, "/", gene_id, "_out.profile"), header = TRUE)

    # Accumulate three variants: raw / z-weighted / sign-weighted
    z_score <- gtex_combine$V1[i]
    keep_gene_expression <- copy(gene_data)
    keep_gene_expression[, SCORE := SCORE * as.numeric(z_score)]
    sign_gene_expression <- copy(gene_data)
    sign_gene_expression[, SCORE := SCORE * sign(as.numeric(z_score))]

    print(head(gene_data))
    print(head(keep_gene_expression))
    print(head(sign_gene_expression))

    if (is.null(all_data)) {
        all_data       <- gene_data
        all_data_keep  <- keep_gene_expression
        all_data_sign  <- sign_gene_expression
    } else {
        all_data <- merge(all_data, gene_data[, .(FID, IID, SCORE)],
                          by = c("FID", "IID"), all = TRUE)
        all_data[, SCORE := rowSums(.SD, na.rm = TRUE),
                 .SDcols = grep("SCORE", names(all_data))]
        all_data <- all_data[, .(FID, IID, PHENO, CNT, CNT2, SCORE)]

        all_data_keep <- merge(all_data_keep, keep_gene_expression[, .(FID, IID, SCORE)],
                               by = c("FID", "IID"), all = TRUE)
        all_data_keep[, SCORE := rowSums(.SD, na.rm = TRUE),
                      .SDcols = grep("SCORE", names(all_data_keep))]
        all_data_keep <- all_data_keep[, .(FID, IID, PHENO, CNT, CNT2, SCORE)]

        all_data_sign <- merge(all_data_sign, sign_gene_expression[, .(FID, IID, SCORE)],
                               by = c("FID", "IID"), all = TRUE)
        all_data_sign[, SCORE := rowSums(.SD, na.rm = TRUE),
                      .SDcols = grep("SCORE", names(all_data_sign))]
        all_data_sign <- all_data_sign[, .(FID, IID, PHENO, CNT, CNT2, SCORE)]
    }

    # Bail early if anything went NA -- means a gene merged badly
    if (any(is.na(all_data))) {
        write(paste("[DEBUG] NA values detected in all_data after iteration", i), stdout())
        print(head(all_data))
        break
    }
    if (any(is.na(all_data_keep))) {
        write(paste("[DEBUG] NA values detected in all_data_keep after iteration", i), stdout())
        print(head(all_data_keep))
        break
    }
}

# Save accumulated per-chromosome scores
all_data_fn      <- paste0(output_dir, "/chr", input_chrm, "_", tissue_type, "_all_data.rds")
all_data_keep_fn <- paste0(output_dir, "/chr", input_chrm, "_", tissue_type, "_all_data_keep.rds")
all_data_sign_fn <- paste0(output_dir, "/chr", input_chrm, "_", tissue_type, "_all_data_sign.rds")
saveRDS(all_data,      file = all_data_fn)
saveRDS(all_data_keep, file = all_data_keep_fn)
saveRDS(all_data_sign, file = all_data_sign_fn)
write("[Script Complete: Data saved successfully]", stdout())

print(head(all_data))
print(head(all_data_keep))
print(head(all_data_sign))
