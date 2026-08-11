#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=3G
#SBATCH --time=5-12:00:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name="bactopia_giardia_2025"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/giardia_assembly/%j_%x.out


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Run bactopia pipeline to perform read QC, assembly, and assembly QC of 
#              the downloaded Giardia duodenalis reads. Run separately for Illumina reads
#              and ONT reads
# Last modified: 2025-11-10
# Usage: sbatch 01_assembly_and_qc.sh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#======================================================================================
#                           Script configuration
#======================================================================================

set -euo pipefail

# Load modules
module load StdEnv/2023 nextflow/25.04.6 apptainer/1.4.5

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source configuration file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo -e "Loading configuration file and defining ENV variables!\n"

script_dir="${SLURM_SUBMIT_DIR:-$PWD}"
if [[ -f "${script_dir}/../config_paths.sh" ]]; then
    source "${script_dir}/../config_paths.sh"
else
    echo "ERROR: '${script_dir}/../config_paths.sh' not found. Run sbatch from script folder!"
    exit 1
fi

#======================================================================================
#                           Bactopia commands (Illumina reads)
#======================================================================================

    # if eagle has low space available in scratch, save temp files in another workdir
nextflow run bactopia/bactopia -r v3.2.0 \
    -profile apptainer,slurm_fir \
    --nfconfig "${script_dir}/configs/bactopia.config" \
    --samples "${processed_data}/bactopia_samplesheets/bactopia_samplesheet_2025.csv" \
    --outdir "${project_outdir}/bactopia_giardia_2025" \
    --shovill_assembler spades \
    --cleanup_workdir \
    --max_genome_size 20000000

#======================================================================================
#                           Bactopia run for hybrid assemblies
#======================================================================================

echo "Running bactopia for hybrid (ONT + Illumina) assemblies"
 
nextflow run bactopia/bactopia -r v3.2.0 \
    -profile apptainer,slurm_fir \
    --sample SAMN12610744 \
    --r1 "${raw_reads_dir}/SRR10007607_1.fastq.gz" \
    --r2 "${raw_reads_dir}/SRR10007607_2.fastq.gz" \
    --ont "${raw_reads_dir}/SRR10007609.fastq.gz" \
    --nfconfig "${script_dir}/configs/bactopia.config" \
    --outdir "${project_outdir}/bactopia_giardia_2025" \
    --short_polish \
    --skip_amr \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 20000000

nextflow run bactopia/bactopia -r v3.2.0 \
    -profile apptainer,slurm_fir \
    --sample SAMN12611599 \
    --r1 "${raw_reads_dir}/SRR10007722_1.fastq.gz" \
    --r2 "${raw_reads_dir}/SRR10007722_2.fastq.gz" \
    --ont "${raw_reads_dir}/SRR10007724.fastq.gz" \
    --nfconfig "${script_dir}/configs/bactopia.config" \
    --outdir "${project_outdir}/bactopia_giardia_2025" \
    --short_polish \
    --skip_amr \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 20000000
