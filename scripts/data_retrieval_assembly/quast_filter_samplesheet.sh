#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name="quast_filter_2025"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/%j_%x.out

#########################################################################################################
#                                           Dependencies
#########################################################################################################

# Load modules
module load StdEnv/2023 apptainer/1.4.5

# Determine script directory (SLURM or local) and source config with paths
set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    script_dir="${SLURM_SUBMIT_DIR}"
else
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi

source "${script_dir}/../../config_paths.sh"


#########################################################################################################
#                                           Command(s)
#########################################################################################################

#--------------------------------------------------------------------------------------------------------
# Assembly QC with quast
#--------------------------------------------------------------------------------------------------------

# Collect all fasta files and create outdir
mkdir -p "${project_outdir}/quast_assembly"

# Appraise quality of all assemblies
apptainer exec "${quast_container}" \
    quast.py \
        --output-dir "${project_outdir}/quast_assembly" \
        --threads 8 \
        --eukaryote \
        --fast \
        "${bactopia_results}"/*/main/assembler/*.fna.gz

#--------------------------------------------------------------------------------------------------------
# Filter Quast output
#--------------------------------------------------------------------------------------------------------

# Filter definition (csvtk): 
    # N50 [col18] must be at least 30,000 [Illumina short reads mostly]
    # n_contigs [col14] less than 1,300 
    # total length [col16] of assembly between 9M and 15M bp [Similar to the refernce genome size of 12M]
cat "${project_outdir}/quast_assembly/transposed_report.tsv" |
    csvtk rename2 --tabs --fuzzy-fields --fields "*" --pattern " " --replacement "_" |
    csvtk rename2 --tabs --fuzzy-fields --fields "*" --pattern "#" --replacement "n" |
    csvtk filter --tabs --filter "N50>30000" | 
    csvtk filter --tabs --filter 'n_contigs<1500' |
    csvtk filter2 --tabs --filter '$Total_length > 9000000 && $Total_length < 15000000' | 
    csvtk cut --tabs --fields Assembly |
    csvtk del-header |
    uniq > "${accessions_dir}/quast_filtered_assemblies.txt"
    

#--------------------------------------------------------------------------------------------------------
# Create samplesheet(s) for nf_chewbbaca
#--------------------------------------------------------------------------------------------------------

# Samplesheet with all assemblies -----------------------------------------------------------------------

# Add header
    # Append all fasta files except failed assemblies
echo "sample,contig" > "${processed_data}/nf_chewbbaca_samplesheets/all_samplesheet_2025.csv"
ls "${bactopia_results}"/*/main/assembler/*.fna.gz |
    grep --invert-match "error.fna.gz$" \
    >> "${processed_data}/nf_chewbbaca_samplesheets/all_samplesheet_2025.csv"

# Samplesheet with HQ assemblies by Quast ---------------------------------------------------------------

# Add header
echo "sample,contig" > "${processed_data}/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv"

# Search each accession in the assembly directory against all assembly paths
while IFS= read -r sample_id; do
    [[ -z "${sample_id}" ]] && continue 
    
    # Match: sample_id after path '/' and flanked by delimiter '.' or '_'
    match=$(grep --max-count=1 "/${sample_id}[._]" "${processed_data}/nf_chewbbaca_samplesheets/all_samplesheet_2025.csv")
    
    # Append sampleid and path if found, or provide error about missing file
    if [[ -n "$match" ]]; then
        echo "${sample_id},${match}" >> "${processed_data}/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv"
    else
        echo "WARNING: No path found for ${sample_id}" >&2
    fi

done < "${accessions_dir}/quast_filtered_assemblies.txt"
