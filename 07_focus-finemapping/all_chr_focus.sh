#!/bin/bash

#SBATCH --partition=shared
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=144G
#SBATCH --time=48:00:00
#SBATCH --export=ALL
#SBATCH --mail-type=END,FAIL          # Send email on job completion or failure
#SBATCH --mail-user=YOUR_EMAIL@example.com  # Your email


# ================================
# Memory & threading guardrails
# ================================

# Prevent BLAS from spawning huge thread pools (CRITICAL)
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Hard cap virtual memory (in KB): 128 GiB
ulimit -v $((128 * 1024 * 1024))




# Ensure a chromosome argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <chromosome>"
    exit 1
fi

CHR=$1  # Assign first command-line argument as chromosome



# Set the path to your sorted_gene.csv file
#sorted_gene_file="${PROJECT_DIR}/FOCUS/21_version/chr${CHR}_genes.csv"
sorted_gene_file="${PROJECT_DIR}/FOCUS/chr${CHR}_genes.csv"


# Initialize the credible set file
#credible_set_file="${PROJECT_DIR}/FOCUS/chr${CHR}_testing_credible_set_pairs.txt"
credible_set_file="${PROJECT_DIR}/FOCUS/chr${CHR}_credible_set_pairs.txt"

> "$credible_set_file"  # Create or empty the file


# Path to the list storing processed gene-tissue pairs
#processed_gene_tissue_list="${PROJECT_DIR}/FOCUS/chr${CHR}_testing_processed_gene_tissue_pairs.txt"
processed_gene_tissue_list="${PROJECT_DIR}/FOCUS/chr${CHR}_processed_gene_tissue_pairs.txt"
> "$processed_gene_tissue_list" # Create or clear the file



# Tracker file to log issues
#tracker_file="${PROJECT_DIR}/FOCUS/chr${CHR}_tracker_testing.txt"
tracker_file="${PROJECT_DIR}/FOCUS/chr${CHR}_tracker.txt"
> "$tracker_file"  # Create or clear the file




# Skip the header and iterate through each line of the CSV file
#tail -n +2 "$sorted_gene_file" | while IFS=$'\t' read -r gene_id V1 p_value PANEL weight WGT ID CHR P0 P1 N SYM START STOP; do
tail -n +2 "$sorted_gene_file" | while IFS=$'\t' read -r gene_id V1 p_value PANEL weight WGT N ID CHR P0 P1 SYM START STOP; do

    
    # Confirm correct parsing
    echo "Processing gene_id: $gene_id, ID: $ID, SYM: $SYM, tissue: $PANEL"
    

    # Define the gene-tissue pair identifier
    gene_tissue_pair="$SYM"
   
    # Check if the gene-tissue pair is already in the credible set file
    if grep -qF "$gene_tissue_pair" "$processed_gene_tissue_list"; then
        echo "Gene-tissue pair ${gene_tissue_pair} already processed. Skipping."
        continue
    fi

    # Define the output .pos file
    #pos_file="${PROJECT_DIR}/FOCUS/output/${SYM}_testing.pos"
    pos_file="${PROJECT_DIR}/FOCUS/output/${SYM}.pos"

    echo "generate the .pos file"	
    # Run the R script to generate the .pos file
    Rscript ${PROJECT_DIR}/FOCUS/Make_Pos_File_forWeights_expanse_v2.R "$gene_id" "$PANEL" "$pos_file"

    # Check if the .pos file exists
    if [[ ! -f "$pos_file" ]]; then
        echo "Missing .pos file for gene_id: ${gene_id}, tissue: ${PANEL}" >> "$tracker_file"
        continue
    fi


    # If the .pos file is successfully created, extract all gene-tissue pairs and add them to the processed list
    if [[ -s "$pos_file" ]]; then
	# Extract gene-tissue pairs from the ID column of the .pos file and append them to the processed list
  	awk 'NR>1 {print $3}' "$pos_file" | while read -r extracted_pair; do
        	# Check if this gene-tissue pair is already in the list before adding
        	if ! grep -qF "$extracted_pair" "$processed_gene_tissue_list"; then
            	    echo "$extracted_pair" >> "$processed_gene_tissue_list"
        	fi
    	done
    fi



    # Change to the ma-focus directory
    cd ${PROJECT_DIR}/FOCUS/ma-focus || { echo "Directory not found"; exit 1; }

    # Initialize Conda for the current shell session
    eval "$(conda shell.bash hook)"

    # Activate the ma-focus environment
    conda activate ma-focus
    
    echo "Create FOCUS DB"
   
 
    # Create FOCUS DB
    #focus import "$pos_file" fusion --name GTEx --assay rnaseq --output "${PROJECT_DIR}/FOCUS/ma-focus/db/fusion_${SYM}_testing_21"
    focus import "$pos_file" fusion --name GTEx --assay rnaseq --output "${PROJECT_DIR}/FOCUS/ma-focus/db/fusion_${SYM}_testing"

    #update db_file name
    #db_file="${PROJECT_DIR}/FOCUS/ma-focus/db/fusion_${SYM}_testing_21.db"
    db_file="${PROJECT_DIR}/FOCUS/ma-focus/db/fusion_${SYM}_testing.db"

   
    # Check if the .db file is empty
    if [ ! -s "$db_file" ]; then
       echo "Skipping: $db_file is empty for gene_id: ${gene_id}, tissue: ${PANEL}" >> "$tracker_file"
       continue
    fi


    
    echo "Finemap"

    CHR_NUM=${CHR#chr}
    # Finemap
    focus finemap ${PROJECT_DIR}/FOCUS/sumstat/chr${CHR_NUM}_formatted_meta_analysis.sumstats.mod2 \
        ${PROJECT_DIR}/TWAS/LDREF/1000G.EUR.${CHR_NUM} \
        "$db_file" \
        --chr "$CHR_NUM" --prior-prob "gencode38" --locations 38:EUR --p-threshold 1 \
        --out "${PROJECT_DIR}/FOCUS/results/${PANEL}/focus_result_${gene_id}"
        #--out "${PROJECT_DIR}/FOCUS/results/${PANEL}/focus_result_${gene_id}_testing"


    
    # Define the path to the FOCUS result file
    #focus_result_file="${PROJECT_DIR}/FOCUS/results/${PANEL}/focus_result_${gene_id}_testing.focus.tsv"
    focus_result_file="${PROJECT_DIR}/FOCUS/results/${PANEL}/focus_result_${gene_id}.focus.tsv"



    # Deactivate the Conda environment
    conda deactivate

    # Return to the original directory
    cd - || exit

    echo "Post-process the finemap results"

    # Post-process the finemap results

    # Check if the FOCUS result file exists and is **not empty**
    if [[ -s "$focus_result_file" ]]; then
        echo "FOCUS result file for ${gene_id}, tissue: ${PANEL} is not empty." >> "$tracker_file"
        # Extract ENS gene ID and tissue where last column == 1
        awk -F'\t' '$NF == 1 {print $4","$5}' "$focus_result_file" | while IFS=, read -r mol_name tissue; do
            # Check for duplicate gene-tissue pairs before adding
            if ! grep -qF "${mol_name},${tissue}" "$credible_set_file"; then
                echo "${mol_name},${tissue}" >> "$credible_set_file"
            fi
        done
    else
        echo "FOCUS result file for ${gene_id}, tissue: ${PANEL} is empty." >> "$tracker_file"
    fi
 




done

