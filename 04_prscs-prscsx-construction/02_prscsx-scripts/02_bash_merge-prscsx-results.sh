#!/bin/bash

# This script merges the PRS-CSx results run on the meta-analysis summary 
# statistics into a single file across all chromosomes for future scoring
# with PLINK; we go through each reference panel's corresponding results 
# separately for the concatenation!

################################################################################
# Note: For these variables in the script, a lot of the paths and names here
#       will be different based on your environment setup! We've removed a lot
#       of the lab-specific path information, but the code logic is the same ...
################################################################################

PROJECT_DIR=/path/to/user/project

INPUTS_DIR=${PROJECT_DIR}/path/to/user/output/from/script/01
OUTPUTS_DIR=${PROJECT_DIR}/path/to/user/output
mkdir -p ${OUTPUTS_DIR}

################################################################################

# - Optional -
# This script contains commands as if we we used the parameter file to 
# evaluate different values of PHI from the previous script! Feel free to
# modify this depending on your needs ... it's just a concatenation :-)

for POPULATION in AFR AMR EAS EUR SAS META; do
    cat "${INPUTS_DIR}"/prscsx_cohort-name_1kg_"${POPULATION}"_pst_eff_a1_b0.5_phi1e-02_chr*.txt \
        > "${OUTPUTS_DIR}"/prscsx_cohort-name_1kg_"${POPULATION}"_pst_eff_a1_b0.5_phi1e-02_allchr.txt
    echo "Done Merging Across Chromosomes for PHI=1e-02 | Population=${POPULATION} | COHORT !"

    cat "${INPUTS_DIR}"/prscsx_cohort-name_1kg_"${POPULATION}"_pst_eff_a1_b0.5_phi1e-04_chr*.txt \
        > "${OUTPUTS_DIR}"/prscsx_cohort-name_1kg_"${POPULATION}"_pst_eff_a1_b0.5_phi1e-04_allchr.txt
    echo "Done Merging Across Chromosomes for PHI=1e-04 | Population=${POPULATION} | COHORT !"

    cat "${INPUTS_DIR}"/prscsx_cohort-name_1kg_"${POPULATION}"_pst_eff_a1_b0.5_phiauto_chr*.txt \
        > "${OUTPUTS_DIR}"/prscsx_cohort-name_1kg_"${POPULATION}"_pst_eff_a1_b0.5_phiauto_allchr.txt
    echo "Done Merging Across Chromosomes for PHI=auto | Population=${POPULATION} | COHORT !"
    
    echo ""
done

echo "Done with script !"
