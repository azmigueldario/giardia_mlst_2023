#!/bin/bash
#SBATCH --mem-per-cpu=5G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=3
#SBATCH --job-name="nfchewbbaca_full"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/giardia_chewbbaca/%j_%x.out

##############################################################################################
#                           Dependencies
##############################################################################################

# set dependencies
module load StdEnv/2023 nextflow/25.10.2 apptainer/1.4.5

# paths
giardia_project_repo="/project/60006/mdprieto/giardia_mlst_2023"
input_samplesheet="${giardia_project_repo}/processed_data/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv"
nf_config_file="${giardia_project_repo}/scripts/nf_chewbbaca/nf_configs/eagle_chewbbaca.config"
outdir="/home/mdprieto/scratch/pilot_results/nf_chewbbaca_pilot_feb2026"

# set manually
nf_chewbbaca_repo="/home/mdprieto/mdp_projects/nf_chewbacca_mlst"
ref_genome_path="/mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/Giardia_GCF000002435_WB_genomic.fna"

##############################################################################################
#                           Nextflow commands
##############################################################################################

rm -f pilot_samplesheet_chewbbaca.csv
# COMMENT OUT FOR FULL RUN
head -n 21 $input_samplesheet > pilot_samplesheet_chewbbaca.csv

if [[ -s pilot_samplesheet_chewbbaca.csv ]]
then
    cat pilot_samplesheet_chewbbaca.csv > active_samplesheet.csv
    rm pilot_samplesheet_chewbbaca.csv
else
    cat $input_samplesheet > active_samplesheet.csv
fi

# nextflow run
nextflow run $nf_chewbbaca_repo/main.nf \
    -resume \
    -profile singularity,slurm \
    -config ${nf_config_file} \
    --input_samplesheet active_samplesheet.csv \
    --outdir ${outdir} \
    --ref_genome ${ref_genome_path} \
    --organism_species "giardia_duodenalis" \
    --eggnog_db "/scratch/group_share/databases/eggnog/eggnog.db" \
    --eggnog_data_dir "/scratch/group_share/databases/eggnog" \
    --eggnog_diamond_db "/scratch/group_share/databases/eggnog/eggnog_proteins.dmnd" \
    --cgMLST_threshold 70

##############################################################################################
#                           Clean-up
##############################################################################################

rm -f active_samplesheet.csv