#!/bin/bash
################################################################################
# Submit one SLURM job per OneK1K cell type to run FUSION.assoc_test.R on the
# meta-analyzed asthma GWAS sumstats. Each job loops chr 1-22 and writes
# per-(cell-type, chromosome) .dat files to TWAS_results_ct/. Existing outputs
# are skipped so re-running fills in only what's missing.
#
# Reads:   TWAS_pos_ct/*.pos  (produced by Make_Pos_File_celltype.R)
# Writes:  TWAS_results_ct/TWAS_metaGWAS_<celltype>_chr<N>.dat
################################################################################

SUMSTATS=${PROJECT_DIR}/TWAS/formatted_meta_analysis.txt
WEIGHTS_DIR=${REF_DIR}/1k1k_cell_type/celltype_weights_v2
LDREF=${PROJECT_DIR}/TWAS/LDREF/1000G.EUR.
SCRIPT=${PROJECT_DIR}/TWAS/fusion_twas/FUSION.assoc_test.R
POSDIR=${PROJECT_DIR}/TWAS/TWAS_pos_ct
OUTDIR=${PROJECT_DIR}/TWAS/TWAS_results_ct
LOGDIR=${PROJECT_DIR}/TWAS/logs_ct

mkdir -p "$OUTDIR" "$LOGDIR"

for posfile in ${POSDIR}/*.pos; do

    celltype=$(basename "$posfile" .pos | sed 's/^twas_//')

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=TWAS_${celltype}
#SBATCH --output=${LOGDIR}/${celltype}.out
#SBATCH --error=${LOGDIR}/${celltype}.err
#SBATCH --partition=shared
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --export=ALL

module load R/4.3.1

echo "Starting ${celltype}"
echo "Node: \$HOSTNAME"

for chr in {1..22}; do
    outfile=${OUTDIR}/TWAS_metaGWAS_${celltype}_chr\${chr}.dat

    if [ -f "\$outfile" ]; then
        echo "Skipping chr \$chr (exists)"
    else
        echo "Running chr \$chr"
        Rscript $SCRIPT \\
            --sumstats $SUMSTATS \\
            --weights $posfile \\
            --weights_dir $WEIGHTS_DIR \\
            --ref_ld_chr $LDREF \\
            --chr \$chr \\
            --out \$outfile
    fi
done

echo "Finished ${celltype}"

EOF

done
