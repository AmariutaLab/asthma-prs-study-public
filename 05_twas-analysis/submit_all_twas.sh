#!/bin/bash
################################################################################
# Submit one SLURM job per GTEx tissue to run FUSION.assoc_test.R on the
# meta-analyzed asthma GWAS sumstats. Each job loops chr 1-22 and writes
# per-(tissue, chromosome) .dat files to TWAS_results/. Existing outputs are
# skipped so re-running fills in only what's missing.
#
# Reads:   TWAS_pos/*.pos  (one row per tissue × gene; produced upstream)
# Writes:  TWAS_results/<tissue>_chr<N>.dat
################################################################################

SUMSTATS=${PROJECT_DIR}/TWAS/formatted_meta_analysis.txt
WEIGHTS_DIR=${REF_DIR}/gtex/weights
LDREF=${PROJECT_DIR}/TWAS/LDREF/1000G.EUR.
SCRIPT=${PROJECT_DIR}/TWAS/fusion_twas/FUSION.assoc_test.R
OUTDIR=${PROJECT_DIR}/TWAS/TWAS_results

mkdir -p logs
mkdir -p $OUTDIR

for posfile in ${PROJECT_DIR}/TWAS/TWAS_pos/*.pos; do

    tissue=$(basename $posfile .pos)

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=TWAS_${tissue}
#SBATCH --output=logs/${tissue}.out
#SBATCH --error=logs/${tissue}.err
#SBATCH --partition=shared
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --export=ALL

module load R/4.3.1

echo "Starting ${tissue}"
echo "Node: \$HOSTNAME"

for chr in {1..22}; do
    outfile=${OUTDIR}/${tissue}_chr\${chr}.dat

    if [ -f "\$outfile" ]; then
        echo "Skipping chr \$chr (exists)"
    else
        echo "Running chr \$chr"
        Rscript $SCRIPT \
            --sumstats $SUMSTATS \
            --weights $posfile \
            --weights_dir $WEIGHTS_DIR \
            --ref_ld_chr $LDREF \
            --chr \$chr \
            --out \$outfile
    fi
done

echo "Finished ${tissue}"

EOF

done

