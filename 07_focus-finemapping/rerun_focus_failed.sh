#!/bin/bash

#SBATCH --partition=shared
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=144G
#SBATCH --time=48:00:00
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@example.com

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

cd ${PROJECT_DIR}/FOCUS/ma-focus

# Ensure conda is active (uncomment if needed, or rely on .bashrc)
# source ~/.bashrc
# conda activate ma-focus

LOGDIR=rerun_logs
mkdir -p $LOGDIR

while read -r line; do
    if [[ "$line" == focus\ finemap* ]]; then
        cmd="$line"
        collecting=1
        args=()
    elif [[ "$collecting" == 1 ]]; then
        if [[ -z "$line" ]]; then
            collecting=0
            # Extract output filename for logging
            out=$(echo "${args[@]}" | sed -n 's/.*--out \([^ ]*\).*/\1/p')
            logfile="$LOGDIR/$(basename $out).log"

            echo "▶ Re-running $out"
            echo "$cmd ${args[*]}" > "$logfile.cmd"

            # Run the command
            $cmd "${args[@]}" &> "$logfile" || echo "FAILED $out" >> failed.list
        else
            # FIX: Remove quotes to split flags and values (e.g., "--chr 22" -> "--chr" "22")
            args+=($line)
        fi
    fi
done < ${PROJECT_DIR}/FOCUS/focus_finemap_before_errors.txt
