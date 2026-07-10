#!/bin/bash
################################################################################
# Single-cell-type FUSION TWAS prototype. Set CT below to the cell type you
# want to score and submit; loops chr 1-22 and writes the per-chr .dat files
# to TWAS_results_ct/. Kept here as a minimal reference for running
# FUSION.assoc_test.R on one CT; the production driver is submit_all_twas_ct.sh
# which iterates every .pos file under TWAS_pos_ct/.
################################################################################

#SBATCH --partition=shared
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --export=ALL



# Load R version 4.3.1
#module load R/4.3.1                  # Adjust this line to load R 4.3.1




CT=b_naive

for chr in {1..22}; do
    Rscript ${PROJECT_DIR}/TWAS/fusion_twas/FUSION.assoc_test.R \
    --sumstats ${PROJECT_DIR}/TWAS/formatted_meta_analysis.txt \
    --weights ${PROJECT_DIR}/TWAS/TWAS_pos_ct/twas_${CT}.pos \
    --weights_dir ${REF_DIR}/1k1k_cell_type/celltype_weights_v2 \
    --ref_ld_chr ${PROJECT_DIR}/TWAS/LDREF/1000G.EUR. \
    --chr $chr \
    --out ${PROJECT_DIR}/TWAS/TWAS_results_ct/TWAS_metaGWAS_${CT}_chr${chr}.dat
done

