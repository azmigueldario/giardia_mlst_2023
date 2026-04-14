#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=3G
#SBATCH --time=5-12:00:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name="bactopia_giardia_2025"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/giardia_assembly/%j_%x.out

#========================================================================================================
#
#                   ASSEMBLY AND QC WITH BACTOPIA 3.0
#
#========================================================================================================

#--------------------------------------------------------------------------------------------------------
#                           Dependencies
#--------------------------------------------------------------------------------------------------------

# Load modules
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
#                           Bactopia commands
#--------------------------------------------------------------------------------------------------------

    # if eagle has low space available in scratch, save temp files in another workdir
nextflow run bactopia/bactopia -r v3.2.0 \
    -profile apptainer,slurm_fir \
    --nfconfig "${script_dir}/bactopia.config" \
    --samples "${processed_data}/bactopia_samplesheets/bactopia_samplesheet_2025.csv" \
    --outdir "${project_outdir}/bactopia_giardia_2025" \
    --shovill_assembler spades \
    --cleanup_workdir \
    --max_genome_size 20000000

