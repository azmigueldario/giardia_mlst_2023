#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                    
#SBATCH --time=1-20:30:00
#SBATCH --account=def-whsiao-ab                 
#SBATCH --cpus-per-task=8
#SBATCH --job-name="sra_download_giardia"           
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=logs_jobs/giardia_assembly/%x_%j.out
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --mail-type=ALL

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
    
    # Prefetch first
    echo "Prefetch of ${SRR}..."
    apptainer exec "${sra_tools_container}" \
        prefetch "${SRR}" \
        --location aws \
        --max-size 100G \
        --output-directory "${scratch_tmp_folder}/sratools_cache" || { 
            echo "ERROR: Failed to prefetch ${SRR}. Skipping..."
            continue 
        }

    # Fasterq-dump
    echo "Fasterq-dump of ${SRR}"
    apptainer exec "${sra_tools_container}" \
        fasterq-dump \
        --split-files \
        --outdir "${raw_reads_dir}" \
        --threads 7 \
        --mem 32GB \
        --temp "${scratch_tmp_folder}/sratools_cache" \
        --force \
        "${scratch_tmp_folder}/sratools_cache/${SRR}/${SRR}.sra"  < /dev/null || { 
            echo "ERROR: Failed to download ${SRR}. Skipping to next..."
            continue 
        }

    echo "Cleaning intermediate file and compressing ${SRR}..."
    rm -rf "${scratch_tmp_folder}/sratools_cache/${SRR}"
    pigz \
        --force \
        --processes 8 \
        "${raw_reads_dir}/${SRR}"*.fastq
        
done < "${accessions_dir}/all_accessions_2025.csv"

#-------------------------------------------------------------------------------------------------------
#                           Reference genome
#-------------------------------------------------------------------------------------------------------

# Download chromosome level assembly of giardia Duodenalis (Illumina + PacBio)
curl \
    --silent \
    --show-error \
    --location \
    --fail \
    https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_duodenalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz \
    > "${reference_genome}/GCF_000002435.2_WB_genomic.fna.gz"

mkdir -p "${reference_genome}/assemblage_B"
curl \
    --silent \
    --show-error \
    --location \
    --fail \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/498/735/GCA_000498735.1_ASM49873v1/GCA_000498735.1_ASM49873v1_genomic.fna.gz" \
    > "${reference_genome}/assemblage_B/GCA_000498735.1_genomic.fna.gz"

