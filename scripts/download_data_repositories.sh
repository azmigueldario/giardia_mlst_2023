#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                    
#SBATCH --time=04:30:00                     
#SBATCH --cpus-per-task=4                   
#SBATCH --job-name="sra_download_biorepository_giardia"           
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/giardia_download_genomes_ena.out   

######################################################################################################

# load necessary modules
module load singularity nextflow

# ENV variables
ACC_LIST="/project/60005/mdprieto/giardia_mlst_2023/accessions/SRR_Acc_List.txt"
BIOR_GENOMES="/project/60005/mdprieto/raw_data/giardia/biorepositories"

# start in scratch
cd /scratch/mdprieto
# download data
nextflow run nf-core/fetchngs -r 1.9 \
    --input "$ACC_LIST" \
    --outdir $BIOR_GENOMES \
    -profile singularity \
    --nf_core_pipeline 'taxprofiler' \
    -resume
    
