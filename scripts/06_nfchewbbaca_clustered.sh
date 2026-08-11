#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=def-whsiao-ab
#SBATCH --job-name="nf_chewbbacca_random_cluster"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --output=logs_jobs/chewbbacca/%j_%x.out


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Run nf_chewbbaca pipeline sensitivity analyses with additional
#              samplesheets (subsampling of dataset) to see if a similar distribution
#              of samples in a minimum spanning tree is shown
# Last modified: 2025-11-10
# Usage: sbatch 05_nf_chewbbaca_clustered.sh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#======================================================================================
#                           Script configuration
#======================================================================================

set -euo pipefail


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source configuration file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo -e "Loading configuration file and defining ENV variables!\n"

script_dir="${SLURM_SUBMIT_DIR:-$PWD}"
if [[ -f "${script_dir}/../config_paths.sh" ]]; then
    source "${script_dir}/../config_paths.sh"
else
    echo "ERROR: ${script_dir}/../config_paths.sh not found. Run sbatch from script folder!"
    exit 1
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                   Define specific paths and load dependencies
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


module load StdEnv/2023 apptainer/1.3.5 nextflow/25.10.2

# Script specific paths
random_samplesheet="${nextflow_samplesheets}/giardia_subsample_random_samplesheet.csv"
minmax_samplesheet="${nextflow_samplesheets}/giardia_subsample_minmax_samplesheet.csv"
nf_config_file="${script_dir}/configs/eagle_chewbbacca.config"
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
        --ref_genome "${assemblage_a_fasta}" \
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
        --ref_genome "${assemblage_a_fasta}" \
        --organism_species "giardia_duodenalis" \
        --eggnog_db "${eggnog_db_dir}/eggnog.db" \
        --eggnog_data_dir "${eggnog_db_dir}" \
        --eggnog_diamond_db "${eggnog_db_dir}/eggnog_proteins.dmnd" \
        --cgMLST_threshold 90 \
        --number_splits 10
