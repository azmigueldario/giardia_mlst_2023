#!/bin/bash
#SBATCH --mem-per-cpu=32G
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name="quast_filter_2025"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out

#########################################################################################################
#                                           Dependencies
#########################################################################################################

# load modules
module load StdEnv/2020 gcc/9.3.0 quast/5.2.0 csvtk/0.23.0

# define environment variables for HPC
project_repo="/project/60006/mdprieto/giardia_mlst_2023"
bactopia_results="/scratch/mdprieto/results/bactopia_giardia_dec2025"
outdir="${project_repo}/output/quast_output"

#########################################################################################################
#                                           Command(s)
#########################################################################################################

# --------------------------------------------- QUAST -------------------------------------------

# Collect all fasta files into a list
ASSEMBLIES=$(find ${bactopia_results} -path "*/main/assembler/*" -type f -name "*fna.gz")

# create directory
mkdir -p ${outdir}

# 
quast.py $ASSEMBLIES \
    --output-dir quast_results \
    --threads 8 \
    --gene-finding \
    --eukaryote \
    --no-html

# -------------------------------------- Filter Quast output ---------------------------------------

csvtk filter2 -t -f "$N50 > 50000 && ${# contigs} < 1300 && ${Total length} > 10000000 && ${Total length} < 14000000" \
    "${outdir}/transposed_report.tsv" | \
    csvtk cut -t -f Assembly > \
        filtered_assemblies.txt


# ------------------------------ Create samplesheet for nf_chewbbaca --------------------------------
