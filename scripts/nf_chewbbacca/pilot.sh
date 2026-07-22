#!/bin/bash
#SBATCH --mem=10G
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=3
#SBATCH --account=def-whsiao-ab
#SBATCH --job-name="nf_chewbbacca_pilot_jun24"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --output=logs_jobs/chewbbacca/%j_%x.out

#----------------------------------------------------------------------------------------
#                           Dependencies and set up
#----------------------------------------------------------------------------------------

# requires apptainer and nextflow > 23
module load StdEnv/2023 apptainer/1.3.5 nextflow/25.10.2

# Import major paths from config files
script_dir="/scratch/mdprieto/repositories/giardia_mlst_2023/scripts/nf_chewbbacca"
source "${script_dir}/../../config_paths.sh"

# Script specific paths
hq_input_samplesheet="${processed_data}/nf_chewbbacca_samplesheets/hq_samplesheet_2025.csv"
nf_config_file="${script_dir}/nf_configs/eagle_chewbbacca.config"

# Define output directory
outdir="${analysis_dir}/nf_chewbbacca/pilot_2026"
mkdir -p "${outdir}"

#----------------------------------------------------------------------------------------
#                           Nextflow commands
#----------------------------------------------------------------------------------------

tmp_samplesheet="/scratch/mdprieto/tmp_samplesheet.csv"
head -n 10 "${hq_input_samplesheet}" > "$tmp_samplesheet"

cd /scratch/mdprieto/ &&
    nextflow run "${nf_chewbbacca_pipeline}/main.nf" \
        -resume \
        -profile apptainer,slurm_fir \
        -w "${nf_work_dir}" \
        -config "${nf_config_file}" \
        --input_samplesheet "${tmp_samplesheet}" \
        --outdir "${outdir}" \
        --ref_genome "${assemblage_A_fasta}" \
        --organism_species "giardia_duodenalis" \
        --eggnog_db "${eggnog_db_dir}/eggnog.db" \
        --eggnog_data_dir "${eggnog_db_dir}" \
        --eggnog_diamond_db "${eggnog_db_dir}/eggnog_proteins.dmnd" \
        --cgMLST_threshold 70 \
        --number_splits 3 \
        --publish_set_data

