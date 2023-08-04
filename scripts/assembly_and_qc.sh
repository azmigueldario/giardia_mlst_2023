#!/bin/bash                                 
#SBATCH --account=rrg-whsiao-ab
#SBATCH --mem-per-cpu=5G                   
#SBATCH --time=12:55:00                     
#SBATCH --cpus-per-task=20                
#SBATCH --job-name="giardia_assembly_and_qc"     
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/%x_%j.out  

###############################################################################################

# load modules
module load nextflow apptainer

# define environment variables
SAMPLESHEET="/project/6056895/mdprieto/giardia_mlst_2023/processed_data/bactopia_samplesheet.csv"
SINGULARITY_LOCAL_CACHE="/project/6007413/cidgoh_share/singularity_imgs"

###############################################################################################

nextflow run bactopia/bactopia -r v2.2.0 \
    -profile singularity \
    -resume \
    --samples $SAMPLESHEET \
    --outdir results/bactopia_giardia \
    --singularity_cache $SINGULARITY_LOCAL_CACHE \
    --shovill_assembler spades \
    --skip_amr 

