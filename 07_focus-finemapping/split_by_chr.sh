#!/bin/bash

#SBATCH --partition=shared
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=1:00:00
#SBATCH --export=ALL



#input_file="${PROJECT_DIR}/FOCUS/sorted_gene_all39Tissue.csv"
#input_file="${PROJECT_DIR}/FOCUS/sorted_gene_all39Tissue_new.csv"
input_file="${PROJECT_DIR}/FOCUS/sorted_gene_all17CT_new.csv"

outdir="${PROJECT_DIR}/FOCUS/chr_genes_17CT"
mkdir -p "$outdir"

# Extract header
header=$(head -n 1 "$input_file")

# Loop over unique chromosomes with 'chr' removed
tail -n +2 "$input_file" | awk -F'\t' '{print $9}' | sort | uniq | while read chr; do
    # Create output file for each chromosome
    echo "$header" > "${outdir}/${chr}_genes.csv"
    awk -F'\t' -v chr="$chr" '$9 == chr' "$input_file" >> "${outdir}/${chr}_genes.csv"
done
