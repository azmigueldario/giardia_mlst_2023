#!/bin/bash
#SBATCH --account=def-whsiao-ab
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name="filtering_clustering_samplesheet"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=logs_jobs/%j_%x.out


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Select HQ bactopia assemblies, sketch genomes with sourmash,
#              apply clustering algorithm using 'scikit', and produce ready-to-run
#              samplesheets for nextflow pipeline            
# Last modified: 2025-11-10
# Usage: sbatch 03_quast_sourmash_tsne.sh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ======================================================================================
#                                           Script configuration
# ======================================================================================

set -euo pipefail

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Source configuration file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo -e "Loading configuration file and defining ENV variables!\n"

script_dir="${SLURM_SUBMIT_DIR:-$PWD}"
if [[ -f "${script_dir}/../config_paths.sh" ]]; then
    source "${script_dir}/../config_paths.sh"
else
    echo "ERROR: ${script_dir}/../config_paths.sh not found. Run sbatch from script folder!"
    exit 1
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load necessary software and define custom local outfiles
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

module load StdEnv/2023 apptainer/1.4.5
all_assemblies_samplesheet="${processed_data}/nf_chewbbaca_samplesheets/all_samplesheet_2025.csv"
hq_assemblies_samplesheet="${processed_data}/nf_chewbbaca_samplesheets/hq_samplesheet_2025.csv"
mkdir -p "$(dirname "${all_assemblies_samplesheet}")"

#======================================================================================
#                           QUAST command(s)
#======================================================================================

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assembly QC with quast
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Filter Quast output
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Filter definition (csvtk): 
#       N50 [col18] must be at least 30,000 [Illumina short reads mostly]
#       n_contigs [col14] less than 1,300 
#       total length [col16] of assembly between 9M and 15M bp [Similar to the refernce genome size of 12M]

cat "${project_outdir}/quast_assembly/transposed_report.tsv" |
    csvtk rename2 --tabs --fuzzy-fields --fields "*" --pattern " " --replacement "_" |
    csvtk rename2 --tabs --fuzzy-fields --fields "*" --pattern "#" --replacement "n" |
    csvtk filter --tabs --filter "N50>30000" | 
    csvtk filter --tabs --filter 'n_contigs<1500' |
    csvtk filter2 --tabs --filter '$Total_length > 9000000 && $Total_length < 15000000' | 
    csvtk cut --tabs --fields Assembly |
    csvtk del-header |
    uniq > "${accessions_dir}/quast_filtered_assemblies.txt"
    

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create samplesheet(s) for nf_chewbbaca with all successfully completed assemblies
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo "sample,contig" > "${all_assemblies_samplesheet}"
for assembly in "${bactopia_results}"/*/main/assembler/*.fna.gz
do
    # Only apply to completed assemblies
    if [[ ! "$assembly" =~ error\.fna\.gz$ ]]
    then   
        sample_name=$(basename "$assembly" .fna.gz)
        echo "${sample_name},${assembly}" \
        >> "${all_assemblies_samplesheet}"
    fi
done

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create samplesheet(s) for nf_chewbbaca with assemblies passing QC filters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

echo "sample,contig" > "${hq_assemblies_samplesheet}"

# Search each accession in the assembly directory against assembly paths
while IFS= read -r sample_id; do
    [[ -z "${sample_id}" ]] && continue 
    
    # Match: sample_id after path '/' and flanked by delimiter '.' or '_'
    match=$(grep --max-count=1 "/${sample_id}[._]" "${all_assemblies_samplesheet}" || true)
    
    # Append sampleid and path if found, or provide error about missing file
    if [[ -n "$match" ]]; then
        echo "${match}" >> "${hq_assemblies_samplesheet}"
    else
        echo "WARNING: No path found for ${sample_id}" >&2
    fi

done < "${accessions_dir}/quast_filtered_assemblies.txt"

#======================================================================================
#                                   Sourmash commands
#======================================================================================

# Create output directory 
sourmash_signatures_dir="${project_outdir}/sourmash/signatures"
sourmash_clean_dir="${project_outdir}/sourmash/signatures_clean"
mkdir -p "${sourmash_signatures_dir}" "${sourmash_clean_dir}"

# Make temporary file for input paths
temp_list="$(mktemp)"

# Parse path column (2) and remove header from list of HQ assemblies
cut \
    --fields=2 \
    --delimiter=, \
    "${hq_assemblies_samplesheet}" |
grep "fna.gz$" \
    > "${temp_list}"

# Sourmash dna sketch scaled 1/1000; specific seed for reproducibility
apptainer exec "${sourmash_container}" \
    sourmash \
        sketch dna \
            --from-file "${temp_list}" \
            --param-string k=51,scaled=1000,seed=1113 \
            --output-dir "${sourmash_signatures_dir}"

# Rename signature identifiers
for signature_file in "${sourmash_signatures_dir}/"*.fna.gz.sig
do
    new_id=$(basename --suffix=".fna.gz.sig" "${signature_file}")
    apptainer exec "${sourmash_container}" \
        sourmash signature rename \
            "${signature_file}" \
            "${new_id}" \
            -o "${sourmash_clean_dir}/${new_id}.fna.gz.sig"
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
            --csv "${sourmash_clean_dir}/../sourmash_distances_hq_genomes.csv" \
            --labels-to "${sourmash_clean_dir}/../sourmash_labels.csv" \
            "${sourmash_clean_dir}/"*.sig
    
rm -r "${temp_list}" "${sourmash_signatures_dir}"

#======================================================================================
#                               t-SNE and HDBSCAN commands
#======================================================================================

# The hyperparameters for t-SNE and HDBSCAN have been previously optimized using the jupyter notebook 
# provided in the repository subfolder ./scripts/clustering_tsne_hdbscan. 
# Here we run t-SNE and HDBSCAN clustering with the selected parameters.

python "${script_dir}/bin/tsne_hdbscan_automated.py" \
    --outfile "${processed_data}/clustered_subsample_list.txt" \
    --random_seed 112233 \
    --min_cluster_size 15 \
    --min_samples 15 \
    "${sourmash_clean_dir}/../sourmash_distances_hq_genomes.csv"