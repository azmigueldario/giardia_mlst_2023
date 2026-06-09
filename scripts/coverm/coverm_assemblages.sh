#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=5G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=10
#SBATCH --job-name="coverm_giardia_assemblages"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/%j_%x.out

#--------------------------------------------------------------------------------------------------------
#                           Dependencies
#--------------------------------------------------------------------------------------------------------

# Load modules
module load StdEnv/2023 apptainer/1.4.5 bwa-mem2/2.2.1

set -euo pipefail

# Determine script directory (SLURM or local) and source config with paths
script_dir="${SLURM_SUBMIT_DIR}"
if [[ ! -f "${script_dir}/../../config_paths.sh" ]]; then
    echo "ERROR: You must cd into the script's directory before running sbatch!"
    exit 1
fi

source "${script_dir}/../../config_paths.sh"

#--------------------------------------------------------------------------------------------------------
#                           Competitive mapping against reference genomes
#--------------------------------------------------------------------------------------------------------

# Create output directory
mkdir -p  "${project_outdir}/coverm"

#  Array with coupled reads (sample1_R1, sample1_R2, sample2_R1, sample2_R2, ... )
coupled_reads=()
for read1 in "${raw_reads_dir}"/*_1.fastq.gz; do
    read2="${read1/_1/_2}"
    coupled_reads+=("$read1" "$read2")
done

# Run alignment for all input reads simultaneously
apptainer exec --cleanenv "${coverm_container}"  \
    coverm genome \
        --genome-fasta-files \
            "${reference_genome_A}" \
            "${reference_genome}/assemblage_B/GCA_000498735.1_ASM49873v1_genomic.fna.gz"  \
        --genome-fasta-extension ".fna.gz" \
        --coupled "${coupled_reads[@]}" \
        --mapper minimap2-sr \
        --methods relative_abundance mean \
        --threads 4 \
        --output-file "${project_outdir}/coverm/competitive_mapping_results.tsv" \
        --min-read-aligned-percent 30
