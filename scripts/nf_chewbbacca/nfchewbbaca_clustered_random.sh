#!/bin/bash
#SBATCH --mem=20G
#SBATCH --time=1-20:00:00
#SBATCH --cpus-per-task=3
#SBATCH --account=def-whsiao-ab
#SBATCH --job-name="nf_chewbbacca_random_cluster"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --output=logs_jobs/chewbbacca/%j_%x.out

#========================================================================================
#
# Core genome MLST pipeline with cross-validation (random sampling of clustered assemblies)
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
input_samplesheet="${processed_data}/nf_chewbbacca_samplesheets/giardia_subsample_random_samplesheet.csv"
nf_config_file="${script_dir}/nf_configs/eagle_chewbbacca.config"

# Define output directory
outdir="${analysis_dir}/nf_chewbbacca/cluster_random_2026"
mkdir -p "${outdir}"

#----------------------------------------------------------------------------------------
#                           Nextflow commands
#----------------------------------------------------------------------------------------

cd /scratch/mdprieto/ &&
    nextflow run "${nf_chewbbacca_pipeline}/main.nf" \
        -resume \
        -profile apptainer,slurm_fir \
        -config "${nf_config_file}" \
        --input_samplesheet "${input_samplesheet}" \
        --outdir "${outdir}" \
        --ref_genome ${reference_genome_A} \
        --organism_species "giardia_duodenalis" \
        --eggnog_db "${eggnog_db_dir}/eggnog.db" \
        --eggnog_data_dir "${eggnog_db_dir}" \
        --eggnog_diamond_db "${eggnog_db_dir}/eggnog_proteins.dmnd" \
        --cgMLST_threshold 90 \
        --number_splits 10
