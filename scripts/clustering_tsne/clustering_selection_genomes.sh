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
#             Clustering and subsampling of BC diversity
##################################################################################################

# With the interactive jupyter notebook, we found a subcluster of BC samples that have 
# more than 97% similarity in the sketchs and that form a distinct subcluster in UMAP.


# As sensitivity analysis, we are producing two samplesheets subsampling (at random and min-max) the 
# diversity in the subcluster

python ${script_dir}/bin/bc_diversity_impact.py \
    --distance_matrix "$DIST_MAT" \
    --metadata "$METADATA" \
    --fasta_dir "$FASTA_DIR" \
    --output_prefix "$OUTPUT_PRE" \
    --fasta_extension .fasta .fa .fna .fasta.gz \
    --agglo_cluster_value 2 \
    --keep 20 \
    --threshold_agglomerative 0.05 \
    --threshold_umap 2.5 \
    --neighbors_umap 15 \
    --seed 1113



    --outfile "${processed_data}/clustered_subsample_list.txt" \
    --random_seed 1113 \
    --min_cluster_size 15 \
    --min_samples 15 \
    ${project_outdir}/sourmash/sourmash_HQ_dist.csv
