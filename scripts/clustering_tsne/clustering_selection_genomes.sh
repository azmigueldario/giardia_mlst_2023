#!/usr/bin/env bash


##################################################################################################
#             Dependencies
##################################################################################################


# load Apptainer software for containers (or have it installed locally)
module load StdEnv/2023 apptainer/1.4.5

# modify project root path as necessary - pointing towards the project repository
project_repo="/project/60006/mdprieto/giardia_mlst_2023"

# manual INPUTs - apptainer image(s) 
sourmash_img="/mnt/cidgoh-object-storage/images/sourmash-4.8.9-hdfd78af_0.img"

# relative paths (taken from previous steps)
repo_bactopia_fastas="${project_repo}/input_data/bactopia_fasta"
quast_filtered_genomes="${project_repo}/processed_data/accessions/quast_filtered_assemblies.txt"
repo_output="${project_repo}/output"


##################################################################################################
#             Sourmash commands
##################################################################################################


mkdir -p ${repo_output}/sourmash/signatures

# 1/1000 scaled dna sketch with a specific seed for reproducibility
apptainer exec ${sourmash_img} \
    sourmash \
        sketch dna \
            $(ls "${repo_bactopia_fastas}"/*.fna.gz | grep -Ef ${quast_filtered_genomes}) \
            --param-string k=51,scaled=1000,seed=1113 \
            --output-dir ${repo_output}/sourmash/signatures

# obtain all vs all distance matrix .csv
apptainer exec ${sourmash_img} \
    sourmash compare \
        --ani \
        --processes 6 \
        --distance-matrix \
        --ksize 51 \
        --dna \
        --csv ${repo_output}/sourmash/sourmash_dist_quast_filtered.csv \
        ${repo_output}/sourmash/signatures/*.fna.gz.sig
    

##################################################################################################
#             t-SNE and HDBSCAN commands
##################################################################################################

# The hyperparameters for t-SNE and HDBSCAN have been previously optimized using the jupyter notebook 
# provided in the repository subfolder ./scripts/clustering_tsne_hdbscan. 
# Here we run t-SNE and HDBSCAN clustering with the selected parameters.

python ${project_root}/scripts/clustering_tsne/bin/tsne_hdbscan_automated.py \
    --outfile "${project_root}/processed_data/clustered_subsample_list.txt" \
    --random_seed 112233 \
    --min_cluster_size 15 \
    --min_samples 15 \
    ${repo_output}/sourmash/sourmash_HQ_dist.csv
