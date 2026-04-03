#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                    
#SBATCH --time=20:30:00                     
#SBATCH --cpus-per-task=8
#SBATCH --job-name="sra_download_giardia"           
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/giardia_assembly/%x_%j.out 

#========================================================================================================
#
#                   DOWNLOAD RAW SEQUENCING READS FROM BIOREPOSITORIES
#
#========================================================================================================

#--------------------------------------------------------------------------------------------------------
#                           Preparation 
#--------------------------------------------------------------------------------------------------------

# Load modules
module load StdEnv/2023 nextflow/25.04.6 apptainer/1.4.5

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
mkdir -p "${scratch_tmp_folder}/sra_tmp"

# Download fastq data and compress it
while read -r SRR; do

    # Download only if not available locally
    if ls "${raw_reads_dir}/${SRR}"*.fastq.gz >/dev/null 2>&1; then
        echo "Skipping ${SRR}: already downloaded and compressed."
        continue
    fi
    
    echo "Downloading ${SRR}..."
    apptainer exec "${sra_tools_container}" \
        fasterq-dump \
        --split-files \
        --outdir "${raw_reads_dir}" \
        --threads 7 \
        --temp "${scratch_tmp_folder}/sra_tmp" \
        --mem 4GB \
        --force \
        "${SRR}" < /dev/null || { 
            echo "ERROR: Failed to download ${SRR}. Skipping to next..."
            continue 
        }

    echo "Compressing ${SRR}..."
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
    > "${reference_genome}/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz"
