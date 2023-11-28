#!/bin/bash
#SBATCH --mem-per-cpu=3G
<<<<<<< HEAD
#SBATCH --time=96:00:00
=======
#SBATCH --time=56:00:00
>>>>>>> 8418fd8 (updates eagle config  with new profile and customized resources)
#SBATCH --cpus-per-task=4
#SBATCH --job-name="giardia_assembly_full"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out

###############################################################################################

# load modules
module load nextflow apptainer

# define environment variables for HPC
SAMPLESHEET="/project/60006/mdprieto/giardia_mlst_2023/processed_data/bactopia_samplesheeet_eagle.csv"
CUSTOM_CONFIG="/project/60006/mdprieto/giardia_mlst_2023/scripts/eagle_bactopia.config "

###############################################################################################

    # eagle has low space available in scratch, save temp files in another temporary workdir
nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    -resume \
    --nfconfig $CUSTOM_CONFIG \
    --samples $SAMPLESHEET \
    --outdir /scratch/mdprieto/results/bactopia_giardia \
    --shovill_assembler spades \
    --skip_amr \
    --long_reads \
    --skip-prokka \
    --skip_mlst \
<<<<<<< HEAD
    --cleanup_workdir 
=======
    --max_memory 90 \
    --cleanup_workdir \
    --slurm_queue 'cpubase_bycore_b1'
>>>>>>> 8418fd8 (updates eagle config  with new profile and customized resources)

