#!/bin/bash
#SBATCH --mem-per-cpu=3G
#SBATCH --time=5-05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name="bactopia_giardia_full"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out

###############################################################################################

# load modules/dependencies
module load nextflow apptainer

# define environment variables for HPC
project_repo="/project/60006/mdprieto/giardia_mlst_2023/"
samplesheet="${project_repo}/processed_data/bactopia_samplesheets/bactopia_samplesheet.csv"
custom_config="${project_repo}scripts/eagle_bactopia.config"
outdir="${HOME}/scratch/results/bactopia_giardia"

###############################################################################################

    # if eagle has low space available in scratch, save temp files in another workdir
nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    -resume \
    --nfconfig ${custom_config} \
    --samples ${samplesheet} \
    --outdir ${outdir} \
    --shovill_assembler spades \
    --skip_amr \
    --long_reads \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 180040666

