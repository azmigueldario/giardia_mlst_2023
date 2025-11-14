#!/usr/bin/env bash


##################################################################################################
#             Dependencies
##################################################################################################


# load Apptainer software for containers (or have it installed locally)
module load apptainer

# modify project root path as necessary - pointing towards the project repository
project_root="/project/60006/mdprieto/giardia_mlst_2023/"

# path to necessary singularity/apptainer image(s)
sourmash_img="/mnt/cidgoh-object-storage/images/sourmash-4.8.9-hdfd78af_0.img"
quast_img="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-quast-5.2.0--py39pl5321heaaa4ec_4.img"

# input paths (taken from previous steps)
bactopia_output_dir="${HOME}/scratch/results/bactopia_giardia"
repo_bactopia_fastas="${project_root}/input_data/bactopia_fasta"

# output paths
repo_output="${project_root}/output"
hq_genomes_list="${project_root}/processed_data/selected_genomes.txt"


##################################################################################################
#             QUAST commands
##################################################################################################


# move all fasta output from bactopia to a single folder in the output directory
mkdir -p ${repo_bactopia_fastas}
find ${bactopia_output_dir} -type f -path "*/assembler/*.fna.gz" -exec cp {} ${repo_bactopia_fastas}/ \;

# run quast for all fasta files (short_reads and hybrid assemblies)
apptainer exec $quast_img quast.py \
    -o ${repo_output}/quast_all_genomes \
    --threads 6 \
    --min-contig 500 \
    --eukaryote \
    --fast \
    $(ls "${repo_bactopia_fastas}"/*.fna)

python ${project_root}/scripts/clustering_tsne/bin/select_genomes_quast.py \
    --outfile ${hq_genomes_list} \
    ${repo_output}/quast_all_genomes/transposed_report.tsv 


##################################################################################################
#             Sourmash commands
##################################################################################################


mkdir -p ${repo_output}/sourmash 

# 1/1000 scaled dna sketch with a specific seed for reproducibility
apptainer exec ${sour_img} sourmash sketch dna \
    $(ls "${repo_bactopia_fastas}"/*.fna | grep -Ef ${hq_genomes_list}) \
    --param-string k=51,scaled=1000,seed=1113 \
    --output-dir ${repo_output}/sourmash 

# obtain all vs all distance matrix .csv
apptainer exec ${sour_img} sourmash compare \
    --ani \
    --processes 6 \
    --distance-matrix \
    --ksize 51 \
    --dna \
    --csv ${repo_output}/sourmash/sourmash_HQ_dist.csv \
    ${repo_output}/sourmash/*.fna.sig
    

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
