#!/bin/bash
#SBATCH --mem-per-cpu=3G
#SBATCH --time=25:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name="pilot_giardia_2025"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/giardia_assembly/%j_%x.out
#SBATCH --mail-user=mprietog@sfu.ca
#SBATCH --mail-type=ALL


#########################################################################################################
#                           Dependencies
#########################################################################################################

# load modules
module load StdEnv/2023 nextflow/25.04.6 apptainer/1.4.5

# define environment variables for HPC
project_repo="/project/60006/mdprieto/giardia_mlst_2023"
samplesheet="${project_repo}/processed_data/bactopia_samplesheets/bactopia_samplesheet_2025.csv"
custom_config="${project_repo}/scripts/data_retrieval_assembly/eagle_bactopia.config"
outdir="${HOME}/scratch/pilot_results/giardia_2025"

#########################################################################################################
#                           Commands
#########################################################################################################

# create pilot samplesheet
head -n 20 $samplesheet > pilot_samplesheet.csv

# if eagle has low space available in scratch, save temp files in another workdir
nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    --nfconfig ${custom_config} \
    --samples pilot_samplesheet.csv \
    --outdir ${outdir} \
    --shovill_assembler spades \
    --cleanup_workdir \
    --max_genome_size 20000000

# delete tmp file
rm pilot_samplesheet.csv