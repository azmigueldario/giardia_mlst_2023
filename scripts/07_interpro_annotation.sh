#!/bin/bash
#SBATCH --mem=64G
#SBATCH --cpus-per-task=12
#SBATCH --time=16:00:00
#SBATCH --job-name="interpro_giardiadb_annotation"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/%j_%x.out

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Verify annotation of loci in final schema
# Last modified: 2026-08-10
# Usage: sbatch interpro_annotation.sh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# =======================================================================================
#                   Configuration
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

module load StdEnv/2023 apptainer/1.3.5 diamond/2.1.22

# Interproscan
nf_chewbbacca_outdir="${analysis_dir}/nf_chewbbacca/full_2026"
input_fasta="${nf_chewbbacca_outdir}/intersect_results/final_results_intersect.fasta"
interpro_outdir="${analysis_dir}/interpro_final_schema"

# GiardiaDB proteome annotation
diamond_outdir="${analysis_dir}/giardiadb_mapping"
diamond_db_reference="${diamond_outdir}/giardiadb_reference"

mkdir -p "${interpro_outdir}/temp"
mkdir -p "${diamond_outdir}"
export APPTAINER_TMPDIR="${interpro_outdir}/temp"

# =======================================================================================
#   Interpro scan command
# =======================================================================================

echo "Running interproscan for fasta input!"

apptainer exec \
    "${interproscan_container}" \
    interproscan.sh \
        --input "${input_fasta}" \
        --output-file-base "${interpro_outdir}/giardia_interpro" \
        --tempdir "${interpro_outdir}/temp" \
        --cpu "${SLURM_CPUS_PER_TASK:-5}" \
        --seqtype n \
        --formats tsv,gff3 \
        --goterms \
        --pathways \
        --disable-precalc \
        --excl-applications SUPERFAMILY,PIRSF,SFLD > "${interpro_outdir}/interproscan_run.log" 2>&1

echo "Finished annotation with interproscan!"

# =======================================================================================
#   Annotate against GiardiaDB proteome
# =======================================================================================

echo "Mapping nucleotides to GiardiaDB reference proteome!"

# Create reference db for search
diamond makedb \
    --in "${giardiadb_proteome}" \
    --db "${diamond_db_reference}"

# 2. Run BLASTx mapping
diamond blastx \
    --db "${diamond_db_reference}" \
    --query "${input_fasta}" \
    --out "${diamond_outdir}/giardiadb_mapping_results.tsv" \
    --threads "${SLURM_CPUS_PER_TASK:-5}" \
    --outfmt 6 qseqid sseqid pident length evalue stitle \
    --max-target-seqs 1 \
    --evalue 1e-5 > "${diamond_outdir}/giardiadb_diamond.log" 2>&1

echo "GiardiaDB mapping complete!"

# =======================================================================================
#   Parse and Merge Annotations (Rscript)
# =======================================================================================

echo "Merging InterProScan and DIAMOND results into final table..."

module purge
module load StdEnv/2023 r/4.5.0

if [[ ! -f "${script_dir}/bin/parse_annotation_tsv.R" ]]; then
    echo "ERROR: R parsing script not found!"
    exit 1
fi

Rscript \
    "${script_dir}/bin/parse_annotation_tsv.R" \
    "${interpro_outdir}/giardia_interpro.tsv" \
    "${diamond_outdir}/giardiadb_mapping_results.tsv" \
    "${analysis_dir}/nfchewbbacca_cleaned_annotations.tsv"