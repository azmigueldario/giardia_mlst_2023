#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=8G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --job-name="coverm_giardia_assemblages"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/%j_%x.out

#--------------------------------------------------------------------------------------------------------
#                           Dependencies
#--------------------------------------------------------------------------------------------------------

module load StdEnv/2023 apptainer/1.4.5 fastp/1.0.1

set -euo pipefail

echo "Determining script directory (SLURM or local) and source 'configuration_file'."

script_dir="${SLURM_SUBMIT_DIR:-$PWD}"
if [[ ! -f "${script_dir}/../../config_paths.sh" ]]; then
    echo "ERROR: You must cd into the script's directory before running sbatch!"
    exit 1
fi
source "${script_dir}/../../config_paths.sh"

#--------------------------------------------------------------------------------------------------------
#                           Local variables
#--------------------------------------------------------------------------------------------------------

export fastp_outdir="${project_outdir}/fastp_clean_reads"
export kraken_outdir="${project_outdir}/kraken2_reports"
export coverm_outdir="${project_outdir}/coverm"

mkdir -p "${fastp_outdir}" "${kraken_outdir}" "${coverm_outdir}"

#--------------------------------------------------------------------------------------------------------
#                           Function for pre-processing and k-mer classification
#--------------------------------------------------------------------------------------------------------

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

    echo "K-mer classification with Kraken2: ${base_name}."
    
    apptainer exec "${kraken2_container}" \
        kraken2 \
            --db "${kraken2_db}" \
            --paired "${clean_r1}" "${clean_r2}" \
            --threads 16 \
            --report "${kraken_outdir}/${base_name}_kraken2_report.txt" \
            --output "${kraken_outdir}/${base_name}_kraken2_output.txt" \
            --gzip-compressed \
            --use-names
done 

echo "Merging all kraken2 reports with krakentools"

kraken_reports=$(ls -1 "${kraken_outdir}"/*_kraken2_report.txt | xargs)
apptainer exec "${krakentools_container}" \
    combine_kreports.py \
        --reports ${kraken_reports} \
        --output "${kraken_outdir}/merged_kraken2_report.txt"

#--------------------------------------------------------------------------------------------------------
#                           Competitive mapping against reference genomes
#--------------------------------------------------------------------------------------------------------

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
            --output-file "${project_outdir}/coverm/competitive_mapping_results.tsv" \
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
            --output-file "${project_outdir}/coverm/competitive_mapping_results_ONT.tsv" \
            --min-read-aligned-percent 10
