#!/bin/bash                                 
#SBATCH --account=rrg-whsiao-ab
#SBATCH --mem-per-cpu=5G                   
#SBATCH --time=02:55:00                     
#SBATCH --cpus-per-task=12                
#SBATCH --job-name="pilot_giardia_assembly_and_qc"     
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/%x_%j.out  

###############################################################################################

# load modules
module load nextflow apptainer

# define environment variables
SAMPLESHEET="/project/6056895/mdprieto/giardia_mlst_2023/processed_data/bactopia_samplesheet.csv"
SINGULARITY_LOCAL_CACHE="/project/6007413/cidgoh_share/singularity_imgs"

###############################################################################################

# bactopia v2.2.0
nextflow run bactopia/main.nf \
    -profile singularity \
    -resume \
    --samples $SAMPLESHEET \
    --outdir results/bactopia_giardia \
    --singularity_cache $SINGULARITY_LOCAL_CACHE \
    --shovill_assembler spades \
    --skip_amr 

