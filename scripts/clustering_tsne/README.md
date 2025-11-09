# Scripts for clustering and selection of representative samples

This step is dependent on the `sci-kit` package for python. The runtime environment can be reproduced by creating a conda environment using the `environment.yml` available in this directory.

## Analysis steps

1. `select_genomes_quast.py` - QC filter to remove highly fragmented and poor quality assemblies from downstream analyses
2. `sourmash_all.sh` - creates a minhash of every genome in the collection to further evaluate distances and perform clustering
2. `tsne_notebook.ipynb` - Interactive Jupyter notebook for embedding and clustering of the population. This notebook produces a simplified samplesheet for the **nf-chewbbaca** pipeline [also available as a python script: `tsne_tuning.py`]