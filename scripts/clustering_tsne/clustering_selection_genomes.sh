#!/usr/bin/env bash


#########################################################################################################
#                                           Dependencies
#########################################################################################################

# Load modules
module load StdEnv/2023 apptainer/1.4.5 python/3.12

# Determine script directory (SLURM or local) and source config with paths
set -euo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    script_dir="${SLURM_SUBMIT_DIR}"
else
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi

source "${script_dir}/../../config_paths.sh"  

##################################################################################################
#             Sourmash commands
##################################################################################################

# Create output directory 
mkdir -p "${project_outdir}/sourmash/signatures"

# Make temporary file for input paths
temp_list="$(mktemp)"

# Parse path column (2) and remove header from list of HQ assemblies
cut \
    --fields=2 \
    --delimiter=, \
     "${processed_data}/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv" |
grep "fna.gz$" \
    > "${temp_list}"


# Sourmash dna sketch scaled 1/1000; specific seed for reproducibility
apptainer exec "${sourmash_container}" \
    sourmash \
        sketch dna \
            --from-file "${temp_list}" \
            --param-string k=51,scaled=1000,seed=1113 \
            --output-dir "${project_outdir}/sourmash/signatures"

# Produce (all vs all) distance matrix csv file
apptainer exec "${sourmash_container}" \
    sourmash \
        compare \
            --ani \
            --processes 6 \
            --distance-matrix \
            --ksize 51 \
            --dna \
            --csv "${project_outdir}/sourmash/sourmash_dist_quast_filtered.csv" \
            "${project_outdir}/sourmash/signatures/"*.fna.gz.sig
    
rm -f "${temp_list}"

##################################################################################################
#             t-SNE and HDBSCAN commands
##################################################################################################

# The hyperparameters for t-SNE and HDBSCAN have been previously optimized using the jupyter notebook 
# provided in the repository subfolder ./scripts/clustering_tsne_hdbscan. 
# Here we run t-SNE and HDBSCAN clustering with the selected parameters.

python ${script_dir}/bin/tsne_hdbscan_automated.py \
    --outfile "${processed_data}/clustered_subsample_list.txt" \
    --random_seed 1113 \
    --min_cluster_size 15 \
    --min_samples 15 \
    ${project_outdir}/sourmash/sourmash_HQ_dist.csv
