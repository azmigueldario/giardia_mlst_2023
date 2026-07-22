#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=2-10:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=def-whsiao-ab
#SBATCH --job-name="nf_chewbbacca_random_cluster"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --output=logs_jobs/chewbbacca/%j_%x.out

#========================================================================================
#
# Core genome MLST pipeline with cross-validation (after clustering)
#
#========================================================================================

#----------------------------------------------------------------------------------------
#                   Source configuration  file
#----------------------------------------------------------------------------------------

set -euo pipefail

script_dir="${SLURM_SUBMIT_DIR}"
if [[ ! -f "${script_dir}/../../config_paths.sh" ]]; then
    echo "ERROR: You must cd into the script's directory before running sbatch!"
    exit 1
fi

source "${script_dir}/../../config_paths.sh"

#----------------------------------------------------------------------------------------
#                   Define specific paths and load dependencies
#----------------------------------------------------------------------------------------

# requires apptainer and nextflow > 23
module load StdEnv/2023 apptainer/1.3.5 nextflow/25.10.2

# Script specific paths
random_samplesheet="${processed_data}/nf_chewbbacca_samplesheets/giardia_subsample_random_samplesheet.csv"
minmax_samplesheet="${processed_data}/nf_chewbbacca_samplesheets/giardia_subsample_minmax_samplesheet.csv"
nf_config_file="${script_dir}/nf_configs/eagle_chewbbacca.config"

# Define output directories
outdir_random="${analysis_dir}/nf_chewbbacca/cluster_random_2026"
outdir_minmax="${analysis_dir}/nf_chewbbacca/cluster_minmax_2026"
mkdir -p "${outdir_random}" "${outdir_minmax}"

#----------------------------------------------------------------------------------------
#                   Nextflow commands - Clustered and subset at random
#----------------------------------------------------------------------------------------

cd /scratch/mdprieto/ &&
    nextflow run "${nf_chewbbacca_pipeline}/main.nf" \
        -resume \
        -profile apptainer,slurm_fir \
        -config "${nf_config_file}" \
        --input_samplesheet "${random_samplesheet}" \
        --outdir "${outdir_random}" \
        --ref_genome "${assemblage_A_fasta}" \
        --organism_species "giardia_duodenalis" \
        --eggnog_db "${eggnog_db_dir}/eggnog.db" \
        --eggnog_data_dir "${eggnog_db_dir}" \
        --eggnog_diamond_db "${eggnog_db_dir}/eggnog_proteins.dmnd" \
        --cgMLST_threshold 90 \
        --number_splits 10


#----------------------------------------------------------------------------------------
#                   Nextflow commands - Clustered and subset to min-max divergence
#----------------------------------------------------------------------------------------

cd /scratch/mdprieto/ &&
    nextflow run "${nf_chewbbacca_pipeline}/main.nf" \
        -resume \
        -profile apptainer,slurm_fir \
        -config "${nf_config_file}" \
        --input_samplesheet "${minmax_samplesheet}" \
        --outdir "${outdir_minmax}" \
        --ref_genome "${reference_genome_A}" \
        --organism_species "giardia_duodenalis" \
        --eggnog_db "${eggnog_db_dir}/eggnog.db" \
        --eggnog_data_dir "${eggnog_db_dir}" \
        --eggnog_diamond_db "${eggnog_db_dir}/eggnog_proteins.dmnd" \
        --cgMLST_threshold 90 \
        --number_splits 10
