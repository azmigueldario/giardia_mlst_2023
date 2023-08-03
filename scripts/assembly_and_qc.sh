#!/bin/bash                                 
#SBATCH --mem-per-cpu=10G                   
#SBATCH --time=06:30:00                     
#SBATCH --cpus-per-task=8                  
#SBATCH --job-name="giardia_assembly_qc_pilot"     
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/%x_%j.out  

###############################################################################################

# load modules
module load nextflow apptainer

# define environment variables
SAMPLESHEET="/project/60006/mdprieto/giardia_mlst_2023/processed_data/bactopia_samplesheet.csv"
SINGULARITY_LOCAL_CACHE="/mnt/cidgoh-object-storage/images"

###############################################################################################

nextflow run bactopia/main.nf \
    -profile singularity \
    -resume \
    --samples $SAMPLESHEET \
    --outdir results/bactopia_giardia \
    --singularity_cache $SINGULARITY_LOCAL_CACHE \
    --shovill_assembler spades \
    --skip_amr 

# -------------- QUAST

quast.py $BC_GENOMES/*.fa \
			-r $GIARDIA_REF/GCF_000002435.2_UU_WB_2.1_genomic.fna \
			-g $GIARDIA_REF/GCF_000002435.2_UU_WB_2.1_genomic.gff \
			-o /scratch/mdprieto/giardia_results/assembly_qc/quast \
			--threads 8

# ------------- CHECKM2