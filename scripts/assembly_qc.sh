#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                   
#SBATCH --time=02:30:00                     
#SBATCH --cpus-per-task=8                  
#SBATCH --job-name="giardia_assembly_qc"     
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/%x_%j.out  

###############################################################################################

# load modules
module load StdEnv/2020  gcc/9.3.0

# define paths
BC_GENOMES="/project/60005/mdprieto/raw_data/giardia/BCCDC/fasta"
GIARDIA_REF="/project/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A"


# -------------- QUAST

quast.py $BC_GENOMES/*.fa \
			-r $GIARDIA_REF/GCF_000002435.2_UU_WB_2.1_genomic.fna \
			-g $GIARDIA_REF/GCF_000002435.2_UU_WB_2.1_genomic.gff \
			-o /scratch/mdprieto/giardia_results/assembly_qc/quast \
			--threads 8

# ------------- CHECKM2