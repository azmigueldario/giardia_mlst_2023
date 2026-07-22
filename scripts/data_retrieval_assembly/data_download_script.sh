#!/bin/bash

#========================================================================================================
#                   DOWNLOAD ILLUMINA READS FROM BIOREPOSITORIES
#========================================================================================================

    # ----------------------------------------------------------------------
    # Dependencies - setup - source paths
    # ----------------------------------------------------------------------

# Load modules
module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9

# Determine script directory
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

# Import configuration
source "${script_dir}/../../config_paths.sh"

# Merge accessions
cat "${accessions_dir}/BCCDC_accessions.txt" "${accessions_dir}/SRR_accessions_2025.csv" > "${accessions_dir}/all_accessions_2025.csv"

# Make directories
mkdir -p "${raw_reads_dir}"
mkdir -p "${reference_genome_dir}"
mkdir -p "${scratch_tmp_folder}/sratools_cache"

    # ----------------------------------------------------------------------
    # Set while loop
    # ----------------------------------------------------------------------

while read -r SRR; do
    # Strip hidden Windows carriage returns
    SRR=$(echo "${SRR}" | tr -d '\r')
    if [ -z "${SRR}" ]; then continue; fi

    if ls "${raw_reads_dir}/${SRR}"*.fastq.gz >/dev/null 2>&1; then
        echo "Skipping ${SRR}: already downloaded."
        continue
    fi

    # ----------------------------------------------------------------------
    # Prefetch with retry
    # ----------------------------------------------------------------------

    max_attempts=2
    attempt=1
    success=false
    echo "Prefetch of ${SRR}..."
    while [ $attempt -le $max_attempts ]; do
        if prefetch "${SRR}" --max-size 100G --output-directory "${scratch_tmp_folder}/sratools_cache" < /dev/null; then
            success=true; break
        else
            echo "WARNING: prefetch attempt $attempt failed for ${SRR}"
            [ $attempt -lt $max_attempts ] && sleep 10m
            attempt=$((attempt + 1))
        fi
    done

    if [ "$success" = false ]; then echo "ERROR: Skipping ${SRR}"; continue; fi

    # ----------------------------------------------------------------------
    # Fasterq-dump of downloaded accessions
    # ----------------------------------------------------------------------

    attempt=1
    success=false
    echo "Fasterq-dump of ${SRR}"
    while [ $attempt -le $max_attempts ]; do
        if fasterq-dump \
                --split-files \
                --outdir "${raw_reads_dir}" \
                --threads 3 \
                --mem 8GB \
                --temp "${scratch_tmp_folder}/sratools_cache" \
                --force "${scratch_tmp_folder}/sratools_cache/${SRR}" < /dev/null; then
            success=true; break
        else
            echo "WARNING: fasterq-dump attempt $attempt failed for ${SRR}."
            [ $attempt -lt $max_attempts ] && sleep 10m
            attempt=$((attempt + 1))
        fi
    done

    if [ "$success" = false ]; then echo "ERROR: Skipping ${SRR}"; continue; fi

    # ----------------------------------------------------------------------
    # Parallelized compression of fastq and cleanup
    # ----------------------------------------------------------------------

    echo "Compressing ${SRR}..."
    rm -rf "${scratch_tmp_folder}/sratools_cache/${SRR}"
    pigz --force --processes 3 "${raw_reads_dir}/${SRR}"*.fastq

done < "${accessions_dir}/all_accessions_2025.csv"

#========================================================================================================
#                   Download reference genome(s)
#========================================================================================================
echo "Downloading reference genome for assemblage A - accession GCF_000002435"
curl \
    --silent \
    --show-error \
    --location \
    --fail \
    https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_duodenalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz \
    > "${assemblage_A_fasta}"

mkdir -p "${reference_genome_dir}/assemblage_B"
echo "Downloading reference genome for assemblage B - accession GCA_000498735"
curl \
    --silent \
    --show-error \
    --location \
    --fail \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/498/735/GCA_000498735.1_ASM49873v1/GCA_000498735.1_ASM49873v1_genomic.fna.gz" \
    > "${reference_genome_dir}/assemblage_B/GCA_000498735.1_genomic.fna.gz"


grep --invert-match --extended-regexp '_[1-2].fastq.gz$'
