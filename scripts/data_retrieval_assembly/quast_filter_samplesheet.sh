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
module load StdEnv/2023 gcc/14.3.0 quast/5.2.0 csvtk/0.23.0 apptainer/1.4.5

# define environment variables for HPC
project_repo="/project/60006/mdprieto/giardia_mlst_2023"
bactopia_results="/scratch/mdprieto/results/giardia_project/bactopia_giardia_feb2025/"
paths="${project_repo}/output/assembly_list.txt"
outdir_quast="${project_repo}/output/quast_output"

# define path to singularity image for quast
quast_container="/scratch/group_share/singularity_imgs/quay.io-quast-5.2.0--py39pl5321h2add14b_1.img"

#########################################################################################################
#                                           Command(s)
#########################################################################################################

# --------------------------------------------- QUAST -------------------------------------------

# Collect all fasta files into a list
ls "${bactopia_results}"/*/main/assembler/*.fna.gz > "$paths"

# create directory
mkdir -p ${outdir}

# to save log in place
cd "$outdir_quast"

# quast for all assemblies
apptainer exec "$quast_container" \
    quast.py $(cat $paths) \
        --output-dir "$outdir_quast" \
        --threads 8 \
        --gene-finding \
        --eukaryote \
        --no-html

# -------------------------------------- Filter Quast output ---------------------------------------

# define output path
filtered_assemblies="${project_repo}/processed_data/accessions/quast_filtered_assemblies.txt"


# csvtk filter definition, columns: N50 is col18, n_contigs is col14, length is col16
csvtk filter2 --tabs \
    --filter '$18 > 50000 && $14 < 1300 && $16 > 9000000 && $16 < 15000000' \
    "${outdir}/transposed_report.tsv" | \
csvtk cut --tabs \
    --fields Assembly > \
    ${filtered_assemblies}


# ------------------------------ Create samplesheet for nf_chewbbaca --------------------------------

# add header to output samplesheet
outfile_samplesheet=""${project_repo}/processed_data/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv""
echo "sample,contig" > "${outfile_samplesheet}"

# loop through each name in the FOFN
while IFS= read -r sample_id; do
    # Skip empty lines
    [[ -z "$sample_id" ]] && continue
    
    # match sample_id after path '/' and must be flanked by delimiter ',' or '_'
    match=$(grep --max-count=1 "/${sample_id}[._]" "$paths")
    
    if [[ -n "$match" ]]; then
        echo "${sample_id},${match}" >> "$outfile_samplesheet"
    else
        echo "WARNING: No path found for ${sample_id}" >&2
    fi

done < "$filtered_assemblies"

# remove file with paths
rm $paths