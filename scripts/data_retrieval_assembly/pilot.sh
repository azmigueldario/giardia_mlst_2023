#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --time=20:30:00
#SBATCH --account=def-whsiao-ab
#SBATCH --cpus-per-task=8
#SBATCH --job-name="sra_download_giardia_pending"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/giardia_assembly/%x_%j.out

#========================================================================================================
#
#                   DOWNLOAD RAW SEQUENCING READS FROM BIOREPOSITORIES
#
#========================================================================================================

#--------------------------------------------------------------------------------------------------------
#                           Preparation
#--------------------------------------------------------------------------------------------------------

# Load modules
module load StdEnv/2023 apptainer/1.4.5

# Determine script directory (SLURM or local)
set -euo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    script_dir="${SLURM_SUBMIT_DIR}"
else
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi

# Import configuration
source "${script_dir}/../../config_paths.sh"

#------------------------------------------------------------------------------------------------------
#                           SRA tools commands
#------------------------------------------------------------------------------------------------------

# Merge all accessions into single file
cat "${accessions_dir}/BCCDC_accessions.txt" "${accessions_dir}/SRR_accessions_2025.csv" > "${accessions_dir}/all_accessions_2025.csv"

# Make directories if necessary
mkdir -p "${raw_reads_dir}"
mkdir -p "${reference_genome}"
mkdir -p "${scratch_tmp_folder}/sratools_cache"

# Download fastq data and compress it
while read -r SRR; do

    # Download only if not available locally
    if ls "${raw_reads_dir}/${SRR}"*.fastq.gz >/dev/null 2>&1; then
        echo "Skipping ${SRR}: already downloaded and compressed."
        continue
    fi

    # ----------------------------------------------------------------------
    # Attempt prefetch with retry
    # ----------------------------------------------------------------------

    max_attempts=2
    attempt=1
    success=false
    echo "Prefetch of ${SRR}..."
    while [ $attempt -le $max_attempts ]; do
        if apptainer exec "${sra_tools_container}" \
                prefetch "${SRR}" \
                --max-size 100G \
                --output-directory "${scratch_tmp_folder}/sratools_cache" < /dev/null; then

            success=true
            break
        else
            echo "WARNING: attempt ${attempt} prefetch for ${SRR} failed"
            if [ $attempt -lt $max_attempts ]; then
                echo "Waiting 10 minutes before retrying..."
                sleep 10m
            fi
            attempt=$((attempt + 1))
        fi
    done

    # Skip if both attempts failed
    if [ "$success" = false ]; then
        echo "ERROR: All prefetch attempts failed for ${SRR}. Skipping to next..."
        continue
    fi

    # ----------------------------------------------------------------------
    # Attempt fasterq-dump with retry
    # ----------------------------------------------------------------------

    echo "Fasterq-dump of ${SRR}"
    max_attempts=2
    attempt=1
    success=false

    # Attempt fasterqdump twice
    while [ $attempt -le $max_attempts ]; do
        if apptainer exec "${sra_tools_container}" \
            fasterq-dump \
            --split-files \
            --outdir "${raw_reads_dir}" \
            --threads 7 \
            --mem 32GB \
            --temp "${scratch_tmp_folder}/sratools_cache" \
            --force \
            "${scratch_tmp_folder}/sratools_cache/${SRR}" < /dev/null; then

            success=true
            break
        else
            echo "WARNING: fasterq-dump attempt $attempt failed for ${SRR}."
            if [ $attempt -lt $max_attempts ]; then
                echo "Waiting 10 minutes before retrying..."
                sleep 10m
            fi
            attempt=$((attempt + 1))
        fi
    done

    # Skip if attempts fail
    if [ "$success" = false ]; then
        echo "ERROR: All fasterq-dump attempts failed for ${SRR}. Skipping to next..."
        continue
    fi

    # ----------------------------------------------------------------------
    # Clean intermediate files and compress raw reads
    # ----------------------------------------------------------------------

    echo "Cleaning intermediate file and compressing ${SRR}..."
    rm -rf "${scratch_tmp_folder}/sratools_cache/${SRR}"
    pigz \
        --force \
        --processes 8 \
        "${raw_reads_dir}/${SRR}"*.fastq

done < "${accessions_dir}/all_accessions_2025.csv"
