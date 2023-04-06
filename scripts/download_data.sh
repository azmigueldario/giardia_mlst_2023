#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                    
#SBATCH --time=04:30:00                     
#SBATCH --cpus-per-task=4                   
#SBATCH --job-name="sra_download"           
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=cfseed_download_data.out   

######################################################################################################

# load necessary modules
module load singularity nextflow

# ENV variables
ACC_LIST="/project/60005/mdprieto/giardia_mlst_2023/accessions/SRR_Acc_List.txt"
OUTPUT_DIR="/project/60005/mdprieto/raw_data/giardia"

# download NCFB
nextflow run nf-core/fetchngs -r 1.9 \
    --input "$ACC_LIST" \
    -profile singularity \
    -resume \
    --outdir $OUTPUT_DIR
