#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=def-whsiao-ab
#SBATCH --job-name="nf_chewbbacca_full"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --output=logs_jobs/chewbbacca/%j_%x.out


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Run nf_chewbbaca pipeline with the selected HQ genomes and performing
#              10 cross-validation steps (loci definition, allele call, 
#              cgMLST definition) before defining the final representative schema#
# Last modified: 2025-11-10
# Usage: sbatch 04_nf_chewbbaca_full.sh
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


# Script specific input and output paths
hq_input_samplesheet="${nextflow_samplesheets}/hq_samplesheet_2025.csv"
nf_config_file="${script_dir}/configs/eagle_chewbbacca.config"
outdir="${analysis_dir}/nf_chewbbacca/full_2026"
mkdir -p "${outdir}"

#======================================================================================
#                          Nextflow commands
#======================================================================================

cd /scratch/mdprieto/ &&
    nextflow run "${nf_chewbbacca_pipeline}/main.nf" \
        -resume \
        -profile apptainer,slurm_fir \
        -config "${nf_config_file}" \
        --input_samplesheet "${hq_input_samplesheet}" \
        --outdir "${outdir}" \
        --ref_genome ${assemblage_a_fasta} \
        --organism_species "giardia_duodenalis" \
        --eggnog_db "${eggnog_db_dir}/eggnog.db" \
        --eggnog_data_dir "${eggnog_db_dir}" \
        --eggnog_diamond_db "${eggnog_db_dir}/eggnog_proteins.dmnd" \
        --cgMLST_threshold 90 \
        --number_splits 10
