#!/bin/bash
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:01:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name="sourmash_all_vs_all"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out

###############################################################################################

# load modules
module load apptainer

# define environment variables for HPC
SOUR_IMG="/mnt/cidgoh-object-storage/images/sourmash-4.8.9-hdfd78af_0.img"
INPUT_DIR="/project/60006/mdprieto/raw_data/bactopia_giardia_fasta"
OUTDIR="../output"
HQ_GENOMES_LIST="../processed_data/selected_genomes.txt"

###############################################################################################

mkdir -p ${OUTDIR}/smash 

# 1/1000 scaled dna sketch with a specific seed
singularity exec $SOUR_IMG sourmash sketch dna \
    $(ls "$INPUT_DIR"/*.fna | grep -Ef $HQ_GENOMES_LIST) \
    --param-string k=51,scaled=1000,seed=1113 \
    --output-dir ${OUTDIR}/smash 

# obtain all vs all distance matrix .csv
singularity exec $SOUR_IMG sourmash compare \
    --ani \
    --processes 6 \
    --distance-matrix \
    --ksize 51 \
    --dna \
    --csv ${OUTDIR}/sourmash_HQ_dist.csv \
    ${OUTDIR}/smash/*.fna.sig