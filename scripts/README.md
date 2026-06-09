# Scripts directory 

# Step 1 - Download and pre-process dataset

1. `./data_retrieval_assembly/download_data_repositories.sh`: Uses nf-core fetchngs and curl to download reference genome and all available NGS data in the INSDC for Giardia duodenalis assemblages A and B

2. `./data_retrieval_assembly/assembly_and_qc.sh`:  Runs bactopia/bactopia 3.0 pipeline for all accessions with short-read data. As _Giardia spp._ genome is mostly redundant across chromosomes, we can use bacterial sequencing assembly tools

3. `./data_retrieval_assembly/hybrid_assembly.sh`: For accessions with short and long sequencing reads, we perform hybrid assembly in Bactopia 3.0



# Step 2 - QC of assembly, clustering, and selecting representative samples

This step is dependent on the `sci-kit` package for python. 

    - The runtime environment can be reproduced by creating a conda environment using the dependencies listed in `./clustering_tsne/tsne_environment.yml`
    - Alternatively, a python `virtual environment` for this step can be created with the dependencies in `./clustering_tsne/tsne_requirements.txt`
    - The jupyter notebook contained in this folder (`./clustering_tsne/tsne_notebook.ipynb`) can be helpful to select the optimal hyperparameters for HDBSCAN clustering.

1. The master script `./clustering_tsne/clustering_selection_genomes.sh` performs all the processes in this directory:

    - Runs `QUAST` for all assemblies resulting from Bactopia
    - Selects high quality genomes based on user-defined parameters
    - Analyzes mash distances among the dataset using `sourmash`
    - Performs clustering and subsampling with `t-SNE` and `HDBSCAN`
    - Creates samplesheets using random subsampling and min-max divergence with user-defined parameters

# Step 3 - nextflow chewbbacca pipeline

These scripts perform gene-by-gene schema definition for the original samplesheet with all assemblies that passed QC, and for the sensitivity analysis using a subset of all available assemblies

1. `./nf_chewbbacca/nfchewbbaca_full.sh`: runs the pipeline for all assemblies that passed QC
2. `./nf_chewbbacca/nfchewbbaca_clustered.sh`: runs the pipeline for a subset of assemblies selected at random from clusters with high genomic similarity. Then it runs the same pipeline with assemblies selected after clustering to maximize divergence in the data.

# Step 4 - Verify alignment to reference assemblages (coverM)
