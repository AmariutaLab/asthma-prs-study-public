#!/usr/bin/env Rscript

# Load the necessary libraries
library(data.table)
# Load necessary library
library(glmnet)
library(dplyr)



# Function to generate the .pos file for a given gene
generate_pos_file <- function(gene_id, tissue, pos_file) {
  # Load the sorted_gene data
  sorted_gene <- fread("${PROJECT_DIR}/FOCUS/sorted_gene_all17CT_new.csv", header = T)
  sorted_gene[, gene_id := sub("\\..*", "", gene_id)]
  
  # Define the focal gene (the one with the lowest p-value)
  focal_gene_id <- strsplit(gene_id, split= "[.]")[[1]][1]
    
  
  # create pos_df
  pos_df <- sorted_gene %>% select("PANEL","WGT","SYM","CHR","P0","P1","N")
  colnames(pos_df)[colnames(pos_df) == "SYM"] <- "ID"  
  
  # remove na values
  pos_df <- na.omit(pos_df)
  
  
  focal_gene_info <- sorted_gene[sorted_gene$gene_id == focal_gene_id | sorted_gene$ID == focal_gene_id, ]
  # Check if focal_gene_info has more than one row
  if (nrow(focal_gene_info) > 1) {
    # Select the row where PANEL matches the tissue
    focal_gene_info <- focal_gene_info[focal_gene_info$PANEL == tissue, , drop = FALSE]
  }
  
  # Get the corresponding row in pos_df
  # Extract chromosome and positions
  focal_chr <- focal_gene_info$CHR
  
  
  focal_start <- as.numeric(focal_gene_info$START) - 1e6
  focal_end <- as.numeric(focal_gene_info$STOP) + 1e6
  
  # Convert to a data frame
  pos_df <- as.data.frame(pos_df)
  
  # Select all genes within ±1 Mb of the focal gene in the same chromosome
  region_genes <- pos_df[
    pos_df$CHR == focal_chr &
      as.numeric(pos_df$P0) <= focal_end &
      as.numeric(pos_df$P1) >= focal_start,
  ]

  # Strip "chr"/"CHR" prefix from all values in the CHR column
  region_genes$CHR <- gsub("^chr", "", tolower(region_genes$CHR))  
  
  # Write the region genes to the .pos file
  fwrite(region_genes, file = pos_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
}





# Main execution: parse command-line arguments and call the function
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: generate_pos_file.R <gene_id> <tissue> <pos_file>")
}

gene_id <- args[1]
tissue <- args[2]
pos_file <- args[3]

generate_pos_file(gene_id, tissue, pos_file)
