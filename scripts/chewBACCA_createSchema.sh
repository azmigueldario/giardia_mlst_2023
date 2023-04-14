#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                   
#SBATCH --time=01:30:00                     
#SBATCH --cpus-per-task=8                  
#SBATCH --job-name="chewBACCA_createSchema_giardia"     
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/%x_%j.out  

###############################################################################################

module load singularity

CHEWBACCA_IMG="/project/cidgoh-object-storage/images/chewbacca_3.1.2.sif"
BC_GENOMES="/project/60005/mdprieto/raw_data/giardia/BCCDC"

export SINGULARITY_BIND="/opt,/scratch,/etc,/project,/mnt"

# run createSchema
singularity exec $CHEWBACCA_IMG chewBBACA.py CreateSchema \
    -i $BC_GENOMES/fasta \
    -o /scratch/mdprieto/giardia_results/chew_createSchema \
    --ptf /project/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/giardia_wb.trn \
    --cpu 8
