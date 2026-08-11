#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=8G
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name="coverm_giardia_assemblages"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/%j_%x.out

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Classify input genomes as either G. duodenalis assemblage A or B. If not
#              enough reads are classified, it will be marked as unknown assemblage. Once
#              the results are parsed, the metadata mastersheet is updated accordingly.
# Last modified: 2026-07-10
# Usage: sbatch 03_quast_sourmash_tsne.sh
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
# Local variables and dependencies
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
module load StdEnv/2023 apptainer/1.4.5 fastp/1.0.1

export coverm_outdir="${project_outdir}/coverm_kraken/coverm"
export fastp_outdir="${project_outdir}/coverm_kraken/fastp_clean_reads"

mkdir -p "${fastp_outdir}" "${kraken_outdir}" "${coverm_outdir}"

#======================================================================================
#                   Function for pre-processing of reads
#======================================================================================

echo "[$(date +'%m-%d %H:%M')] Starting fastp and kraken2..."

for read1 in "${raw_reads_dir}"/*_1.fastq.gz
do
    read2="${read1/_1.fastq.gz/_2.fastq.gz}"
    base_name=$(basename "$read1" _1.fastq.gz)

    clean_r1="${fastp_outdir}/${base_name}_1-clean.fastq.gz"
    clean_r2="${fastp_outdir}/${base_name}_2-clean.fastq.gz"

    echo "Processing with fastp: ${base_name}."

    fastp \
        --in1 "${read1}" --in2 "${read2}" \
        --out1 "${clean_r1}" --out2 "${clean_r2}" \
        --json "${fastp_outdir}/${base_name}_fastp.json" \
        --thread 4 \
        --detect_adapter_for_pe

done

#======================================================================================
#                   Competitive mapping against reference genomes
#======================================================================================

echo "Prepare array with coupled reads: ['sample1_R1', 'sample1_R2', 'sample2_R1', 'sample2_R2', ... ]"

coupled_reads=()
for read1 in "${fastp_outdir}"/*_1-clean.fastq.gz; do
    read2="${read1/_1-clean.fastq.gz/_2-clean.fastq.gz}"
    coupled_reads+=("$read1" "$read2")
done

echo "Run competitive alignment with CoverM for Illumina reads."

apptainer exec --cleanenv "${coverm_container}"  \
    coverm \
        genome \
            --genome-fasta-files \
                "${assemblage_A_fasta}" \
                "${reference_genome_dir}/assemblage_B/GCA_000498735.1_genomic.fna.gz"  \
            --coupled "${coupled_reads[@]}" \
            --mapper minimap2-sr \
            --methods relative_abundance mean covered_bases \
            --threads 16 \
            --output-file "${coverm_outdir}/competitive_mapping_results.tsv" \
            --min-read-aligned-percent 10

echo "Run competitive alignment with CoverM for ONT reads."

apptainer exec --cleanenv "${coverm_container}"  \
    coverm \
        genome \
            --genome-fasta-files \
                "${assemblage_A_fasta}" \
                "${reference_genome_dir}/assemblage_B/GCA_000498735.1_genomic.fna.gz"  \
            --single \
                "${raw_reads_dir}/SRR10007609.fastq.gz" \
                "${raw_reads_dir}/SRR10007724.fastq.gz" \
            --mapper minimap2-ont \
            --methods relative_abundance mean covered_bases \
            --threads 16 \
            --output-file "${coverm_outdir}/competitive_mapping_results_ONT.tsv" \
            --min-read-aligned-percent 10

#======================================================================================
#                       Parse coverM mapping results - update metadata
#======================================================================================

module purge
module load python/3.11
source "${clustering_env}/bin/activate"

python3 "${script_dir}/bin/map_giardia_assemblages.py" \
    --abundance "${coverm_outdir}/competitive_mapping_results_ONT.tsv" \
                "${coverm_outdir}/competitive_mapping_results.tsv" \
    --master "${processed_data}/metadata/merged_SraPrystajecky.csv" \
    --sample-col "run" \
    --output "${processed_data}/metadata/final_SraPrystajecky.csv"
