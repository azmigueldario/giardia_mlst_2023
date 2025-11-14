# Scripts for clustering and selection of representative samples

This step is dependent on the `sci-kit` package for python. The runtime environment can be reproduced by creating a conda environment using the `environment.yml` available in this directory.

The master script `clustering_selection_genomes.sh` extracts all fasta files resulting from the Bactopia run into a single folder, runs `QUAST` and selects high quality genomes based on pre-defined parameters,
analyzes mash distances among the dataset using `sourmash`, and performs clustering and subsampling with `t-SNE` and `HDBSCAN`

The jupyter notebook contained in this folder (`tsne_notebook.ipynb`) can be helpful to select the optimal hyperparameters for HDBSCAN clustering.