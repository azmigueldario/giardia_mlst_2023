#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --job-name="giardia_assembly_qc_full"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%x_%j.out

###############################################################################################

# load modules
module load nextflow apptainer

# define environment variables for HPC
SAMPLESHEET="/project/60006/mdprieto/giardia_mlst_2023/processed_data/eagle_bactopia.csv"
CUSTOM_CONFIG="/project/60006/mdprieto/giardia_mlst_2023/scripts/eagle/eagle.config"

###############################################################################################

nextflow run bactopia/bactopia -r v2.2.0 \
-profile singularity \
-resume \
--nfconfig $CUSTOM_CONFIG \
--samples $SAMPLESHEET \
--outdir results/bactopia_giardia \
--shovill_assembler spades \
--skip_amr

