#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name="filtering_clustering_samplesheet"
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
#                                          QUAST command(s)
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
# Create samplesheet(s) for nf_chewbbaca with all assemblies
#--------------------------------------------------------------------------------------------------------

# Add header
    # Append all fasta files except failed assemblies
echo "sample,contig" > "${processed_data}/nf_chewbbaca_samplesheets/all_samplesheet_2025.csv"
ls "${bactopia_results}"/*/main/assembler/*.fna.gz |
    grep --invert-match "error.fna.gz$" \
    >> "${processed_data}/nf_chewbbaca_samplesheets/all_samplesheet_2025.csv"

#--------------------------------------------------------------------------------------------------------
# Create samplesheet(s) for nf_chewbbaca with assemblies passing QC filters
#--------------------------------------------------------------------------------------------------------

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


##################################################################################################
#                                   Sourmash commands
##################################################################################################

# Create output directory 
mkdir -p "${project_outdir}/sourmash/signatures"
mkdir -p "${project_outdir}/sourmash/signatures_clean"

# Make temporary file for input paths
temp_list="$(mktemp)"

# Parse path column (2) and remove header from list of HQ assemblies
cut \
    --fields=2 \
    --delimiter=, \
     "${processed_data}/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv" |
grep "fna.gz$" \
    > "${temp_list}"


# Sourmash dna sketch scaled 1/1000; specific seed for reproducibility
apptainer exec "${sourmash_container}" \
    sourmash \
        sketch dna \
            --from-file "${temp_list}" \
            --param-string k=51,scaled=1000,seed=1113 \
            --output-dir "${project_outdir}/sourmash/signatures"

# Rename signature identifiers
for signature_file in "${project_outdir}/sourmash/signatures/"*.fna.gz.sig
do
    new_id=$(basename --suffix=".fna.gz.sig" "${signature_file}")
    apptainer exec "${sourmash_container}" \
        sourmash signature rename \
            "${signature_file}" \
            "${new_id}" \
            -o "${project_outdir}/sourmash/signatures_clean/${new_id}.fna.gz.sig"
done

# Produce (all vs all) distance matrix csv file
apptainer exec "${sourmash_container}" \
    sourmash \
        compare \
            --ani \
            --processes 6 \
            --distance-matrix \
            --ksize 51 \
            --dna \
            --csv "${project_outdir}/sourmash/sourmash_distances_hq_genomes.csv" \
            --labels-to "${project_outdir}/sourmash/sourmash_labels.csv" \
            "${project_outdir}/sourmash/signatures_clean/"*.sig
    
rm "${temp_list}" "${project_outdir}/sourmash/signatures"

##################################################################################################
#                               t-SNE and HDBSCAN commands
##################################################################################################

# The hyperparameters for t-SNE and HDBSCAN have been previously optimized using the jupyter notebook 
# provided in the repository subfolder ./scripts/clustering_tsne_hdbscan. 
# Here we run t-SNE and HDBSCAN clustering with the selected parameters.

python ${project_root}/scripts/clustering_tsne/bin/tsne_hdbscan_automated.py \
    --outfile "${project_root}/processed_data/clustered_subsample_list.txt" \
    --random_seed 112233 \
    --min_cluster_size 15 \
    --min_samples 15 \
    ${project_outdir}/sourmash/sourmash_HQ_dist.csv