#!/bin/bash
#SBATCH --mem-per-cpu=3G
#SBATCH --time=1-05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name="bactopia_hybrid_2025"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/giardia_assembly/%j_%x.out

#========================================================================================================
#
#                       HYBRID ASSEMBLY WITH BACTOPIA 3.0
#
#========================================================================================================

#--------------------------------------------------------------------------------------------------------
#                           Dependencies
#--------------------------------------------------------------------------------------------------------

# Load modules/dependencies
module load StdEnv/2023 nextflow/25.04.6 apptainer/1.4.5

# Determine script directory (SLURM or local) and source config with paths
set -euo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    script_dir="${SLURM_SUBMIT_DIR}"
else
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
source "${script_dir}/../../config_paths.sh"


#--------------------------------------------------------------------------------------------------------
#                           Bactopia run for hybrid assemblies
#--------------------------------------------------------------------------------------------------------

nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    --sample SAMN12610744 \
    --r1 "${raw_reads_dir}/SRR10007607_1.fastq.gz" \
    --r2 "${raw_reads_dir}/SRR10007607_2.fastq.gz" \
    --ont "${raw_reads_dir}/SRR10007609.fastq.gz" \
    --nfconfig "${script_dir}/eagle_bactopia.config" \
    --outdir "${project_outdir}/bactopia_giardia_2025" \
    --short_polish \
    --skip_amr \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 180040666

nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    --sample SAMN12611599 \
    --r1 "${raw_reads_dir}/SRR10007722_1.fastq.gz" \
    --r2 "${raw_reads_dir}/SRR10007722_2.fastq.gz" \
    --ont "${raw_reads_dir}/SRR10007724.fastq.gz" \
    --nfconfig "${script_dir}/eagle_bactopia.config" \
    --outdir "${project_outdir}/bactopia_giardia_2025" \
    --short_polish \
    --skip_amr \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 180040666
