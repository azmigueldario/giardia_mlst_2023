#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:30:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name="giardia_assembly_qc_pilot"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%x_%j.out

###############################################################################################

# load modules
module load nextflow apptainer

# define environment variables for HPC (replace as necessary)
SAMPLESHEET="/project/60006/mdprieto/giardia_mlst_2023/processed_data/pilot_bactopia_samplesheeet.csv"
CUSTOM_CONFIG="/project/60006/mdprieto/giardia_mlst_2023/scripts/eagle_bactopia.config "

###############################################################################################

nextflow run bactopia/bactopia -r v3.0.0 \
-profile singularity \
--nfconfig $CUSTOM_CONFIG \
--samples $SAMPLESHEET \
--outdir /scratch/mdprieto/results/bactopia_giardia \
--shovill_assembler spades \
--skip_amr \
--long_reads \
--skip-prokka \
--skip_mlst


