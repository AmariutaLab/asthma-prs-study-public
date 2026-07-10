#!/bin/bash
################################################################################
# Single-tissue FUSION TWAS prototype (whole_blood). Loops chr 1-22 and writes
# the per-chr .dat files into the working directory. Kept here as a minimal
# reference for running FUSION.assoc_test.R interactively / on one tissue;
# the production driver is submit_all_twas.sh which iterates every .pos file
# under TWAS_pos/ and writes into TWAS_results/.
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




for chr in {1..22}; do
    Rscript ${PROJECT_DIR}/TWAS/fusion_twas/FUSION.assoc_test.R \
    --sumstats ${PROJECT_DIR}/TWAS/formatted_meta_analysis.txt \
    --weights ${PROJECT_DIR}/TWAS/TWAS_pos/twas_whole_blood.pos \
    --weights_dir ${REF_DIR}/gtex/weights \
    --ref_ld_chr ${PROJECT_DIR}/TWAS/LDREF/1000G.EUR. \
    --chr $chr \
    --out ${PROJECT_DIR}/TWAS/TWAS_metaGWAS_blood_chr${chr}.dat
done

