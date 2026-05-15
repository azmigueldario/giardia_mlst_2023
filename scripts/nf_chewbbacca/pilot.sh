#!/bin/bash
#SBATCH --mem=80G
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=8
#SBATCH --account=def-whsiao-ab
#SBATCH --job-name="nf_chewbbacca_pilot"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --output=logs_jobs/chewbbacca/%j_%x.out

#========================================================================================
#
# Core genome MLST pipeline with cross-validation
#
#========================================================================================

#----------------------------------------------------------------------------------------
#                   Set dependencies, paths, and containers
#----------------------------------------------------------------------------------------

# requires apptainer and nextflow > 23
module load StdEnv/2023 apptainer/1.3.5 nextflow/25.10.2

# Determine script directory (SLURM or local)
set -euo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    script_dir="${SLURM_SUBMIT_DIR}"
else
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi

# Import major paths from config files
source "${script_dir}/../../config_paths.sh"

# Script specific paths
input_samplesheet="${processed_data}/nf_chewbbacca_samplesheets/hq_samplesheet_2025.csv"
nf_config_file="${script_dir}/nf_configs/eagle_chewbbacca.config"
outdir="${analysis_dir}/nf_chewbbacca/pilot"
mkdir -p "${outdir}"

#----------------------------------------------------------------------------------------
#                           Nextflow commands
#----------------------------------------------------------------------------------------

tmp_samplesheet="${outdir}/pilot_samplesheet.csv"
cat "${input_samplesheet}" | head -n 15 > "${tmp_samplesheet}"

# nextflow run
cd /scratch/mdprieto/ &&
    nextflow run "${nf_chewbbacca_pipeline}/main.nf" \
        -resume \
        -profile apptainer,slurm_fir \
        -config "${nf_config_file}" \
        --input_samplesheet "${tmp_samplesheet}" \
        --outdir "${outdir}" \
        --ref_genome ${reference_genome_A} \
        --organism_species "giardia_duodenalis" \
        --eggnog_db "${eggnog_db_dir}/eggnog.db" \
        --eggnog_data_dir "${eggnog_db_dir}" \
        --eggnog_diamond_db "${eggnog_db_dir}/eggnog_proteins.dmnd" \
        --cgMLST_threshold 70 \
        --number_splits 4
