#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                   
#SBATCH --time=00:30:00                     
#SBATCH --cpus-per-task=4                  
#SBATCH --job-name="prodigal_training_giardia"     
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/prodigal_training.out  

###############################################################################################

# reference genome for Giardia intestinalis assemblage A
GIARDIA_REF="/project/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A"

# singularity prodigal image
PRODIGAL_IMG="/project/cidgoh-object-storage/images/prodigal_2.6.3.sif"

# load singularity to instance
module load singularity
export SINGULARITY_BIND="/opt,/scratch,/etc,/project,/mnt"

# run prodigal training
singularity exec $PRODIGAL_IMG prodigal \
    -i $GIARDIA_REF/GCF_000002435.2_UU_WB_2.1_genomic.fna \
    -t $GIARDIA_REF/giardia_wb.trn \
    -p single 