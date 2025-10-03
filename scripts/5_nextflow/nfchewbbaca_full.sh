#!/bin/bash
#SBATCH --mem-per-cpu=5G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=3
#SBATCH --job-name="nfchewbbaca_full"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out

############################################################

# Requires apptainer and nextflow
module load apptainer
source /project/share/tools/anaconda3/etc/profile.d/conda.sh
conda activate nf-core-tools

############################################################


NF_CHEW="/project/60006/mdprieto/nf_chewbacca_mlst"

cd ~/scratch/ &&
    nextflow run $NF_CHEW/main.nf \
    -resume \
    -profile singularity \
    -c /project/60006/mdprieto/nf_chewbacca_mlst/test/eagle.config \
    --input_samplesheet /project/60006/mdprieto/nf_chewbacca_mlst/test/full_HQ_samplesheet.csv \
    --outdir /scratch/mdprieto/results/nf_chewbbaca/apr2025/full \
    -with-trace