#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name="classic_schemas_analysis"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/%j_%x.out


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Evaluates how classic MLST schemas (using three of six CDS) can classify
#              the full study dataset. First, it downloads all necessary data and an
#              outgroup assembly. Then, it applies a function that employs chewBBACA to
#              map the loci in the schema, create a multiple sequence alignment, and a
#              phylogenetic tree with IQTREE.
# Last modified: 2026-06-10
# Usage: sbatch 04_classic_schema_analysis.sh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# =======================================================================================
#
# Configuration
#
# =======================================================================================

set -euo pipefail

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source configuration file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo -e "Loading configuration file and defining ENV variables!\n"
script_dir="${SLURM_SUBMIT_DIR:-$PWD}"
if [[ -f "${script_dir}/../config_paths.sh" ]]; then
    source "${script_dir}/../config_paths.sh"
else
    echo "ERROR: ${script_dir}/../config_paths.sh not found. Run sbatch from script folder!"
    exit 1
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load dependencies and define ENV variables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

module load StdEnv/2023 apptainer/1.3.5

mlst_fasta_three_loci="${input_dir}/classic_schemas/three_loci"
mlst_fasta_six_loci="${input_dir}/classic_schemas/six_loci"
outdir_classic_three="${analysis_dir}/classic_schemas/three_loci"
outdir_classic_six="${analysis_dir}/classic_schemas/six_loci"
input_filepaths_list="${accessions_dir}/fofn_HQclustered_classicMLST.txt"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run helper script that downloads data for this analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo -e "Running helper script that downloads necessary data!\n"
bash "${script_dir}/bin/classic_schemas_data_preparation.sh"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parse samplesheet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo -e "Preparing samplesheet and decompressing fasta files!\n"

# Add outgroup reference assembly
echo "${input_dir}/classic_schemas/GCA_006247105-1_1_genomic.fna.gz" \
    > "${input_filepaths_list}"

# Select only the fasta files that passed QC
tail --lines=+2 "${processed_data}/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv" \
    | cut --fields=2 --delimiter=, \
    >> "${input_filepaths_list}"

# Decompress fasta files if necessary, fix naming in FOFN
grep '\.gz$' "${input_filepaths_list}" \
    | xargs --no-run-if-empty gunzip --keep \
    || true

sed --in-place 's/\.gz$//' "${input_filepaths_list}"

# =======================================================================================
#
#  Define function:
#   -AlleleCalling, Multiple Sequence Alignment, and phylogenetic tree.
#
# =======================================================================================

run_classic_mlst_pipeline(){

    # Define required inputs
    local schema_dir="$1"
    local outdir="$2"
    local prefix="$3"

    # Create root outdir
    rm -rf "${outdir}"
    mkdir -p "${outdir}"

    # 1. Use the downloaded loci to create an chewBBACA schema
    apptainer exec "${chewbbaca_container}" chewBBACA.py PrepExternalSchema \
        --schema-directory "${schema_dir}" \
        --translation-table 1 \
        --output-directory "${outdir}/adapted_schema" \
        --cpu-cores 6

    # 2. Perform allele calling using the adapted schema
    apptainer exec "${chewbbaca_container}" chewBBACA.py AlleleCall \
        --input-files "${input_filepaths_list}" \
        --schema-directory "${outdir}/adapted_schema" \
        --translation-table 1 \
        --output-directory "${outdir}/AlleleCall" \
        --cpu-cores 6

    # 3. Extract the matching alleles as fasta files
    apptainer exec "${chewbbaca_container}" chewBBACA.py GetAlleles \
        --input-file "${outdir}/AlleleCall/results_alleles.tsv"  \
        --schema-directory "${outdir}/adapted_schema" \
        --translation-table 1 \
        --output-directory "${outdir}/GetAlleles" \
        --cpu-cores 6

    # 4. Perform multiple sequence alignment (MSA) for each locus
    echo -e "\n\nMultiple sequence alignment of Loci from MLST schema.\n\n"
    mkdir -p "${outdir}/multiple_alignment"
    for fasta_file in "${outdir}"/GetAlleles/fastas/*.fasta
    do
        filename=$(basename -- "${fasta_file}" ".fasta")
        apptainer exec "${chewbbaca_container}" mafft \
            --globalpair \
            --maxiterate 1000 \
            --thread $(( ${SLURM_CPUS_PER_TASK:-4} - 1 )) \
            "${fasta_file}" > "${outdir}/multiple_alignment/${filename}_aligned.fasta"
    done

    # 5. Combine the multiple sequence alignments into a single file
    echo -e "\n\nMerging all multiple sequence alignment files into single multifasta.\n\n"
    apptainer exec "${seqkit_container}" \
        seqkit concat \
            --full \
            --fill "-" \
            --id-regexp '(.*(?:SRR|ERR|SAMN|GCA_?)[0-9]+)' \
            --line-width 80 \
            --threads 6 \
            --out-file "${outdir}/multiple_alignment/merged_MSA.fasta" \
            "${outdir}"/multiple_alignment/*_aligned.fasta

    # 6. Create simple phylogenetic tree (IQTREE2) with automated model selection '-m',
    #       1000 ultrafast bootstraps '-B' , and 1000 replicates 'alrt',
    #

    echo -e "\n\nCalculating phylogenetic tree with IQTREE3.\n\n"
    mkdir -p "${outdir}/iqtree"
    apptainer exec "${iqtree_container}" iqtree \
        -s "${outdir}/multiple_alignment/merged_MSA.fasta" \
        -m TEST \
        -T AUTO \
        -B 1000 \
        --pathogen \
        --alrt 1000  \
        --prefix "${outdir}/iqtree/${prefix}" \
        --verbose
}

# =======================================================================================
#
#   Execution
#
# =======================================================================================

echo -e "Running pipeline of analysis for three-loci schema!\n"

run_classic_mlst_pipeline \
    "${mlst_fasta_three_loci}" \
    "${outdir_classic_three}" \
    "three_loci"

echo -e "Running pipeline of analysis for six-loci schema!\n"

run_classic_mlst_pipeline \
    "${mlst_fasta_six_loci}" \
    "${outdir_classic_six}" \
    "six_loci"

echo -e "Cleaning decompressed files if there is gzipped backup available!\n"

while read -r uncompressed_file; do
    if [[ -f "${uncompressed_file}.gz" ]]; then
        rm "${uncompressed_file}"
    fi
done < "${input_filepaths_list}"
