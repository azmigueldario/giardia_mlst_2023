#!/bin/bash
#SBATCH --mem-per-cpu=5G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=3
#SBATCH --job-name="bwa_full"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out

############################################################

# Requires apptainer container functionality 
module load apptainer

    # add conda functionality in compute node for nextflow
source /project/share/tools/anaconda3/etc/profile.d/conda.sh
conda activate nf-core-tools

############################################################

cd ~/scratch/test && \
nextflow run /home/mdprieto/scratch/bwa_align_stats/main.nf \
    -profile singularity,slurm \
    -resume \
    -params-file /scratch/mdprieto/bwa_align_stats/assets/eagle_configs/params_eagle.json \
    -config /scratch/mdprieto/bwa_align_stats/assets/eagle_configs/eagle_cidgoh.config
