# Scripts directory 

## Preparation

The root folder of the **repository** as well as the necessary path(s) to container(s) and dependencies shouuld be modified in the `../config_paths.sh` file

All the scripts are dependent on the configuration file. 

## Step 1 - Download data and pre-process reads

1. `./00_data_download_script.sh`: ideally run in a persistent window (`nohup` or `tmux`), uses the SRA toolkit to download all reads and reference genomes into subfolders of the repository

2. `./01_assembly_and_qc.sh`:  runs bactopia/bactopia 3.0 pipeline for all accessions (separately for Illumina and hybrid assemblies) with short-read data. We apply bacterial assembly tools as _Giardia duodenalis_ genome is mostly redundant across chromosomes.

3. `./02_coverm_assemblages.sh`: classifies the reads as either assemblage A, assemblage B, or unknown (mostly because of limited mapping) using competitive mapping against the reference genomes of the assemblages

## Step 2 - Filtering and parsing assemblies

- This step is dependent on python machine learning tools like `sci-kit`. You can reproduce the **virtual environment** with the dependencies in `./configs/env_recipes/tsne_requirements.txt`. Once the `venv` is created, specify its path in the primary configuration file (`../config_paths.sh`).

    - The jupyter notebook contained in the bin folder [`bin/bc_diversity_impact.ipynb`] helps select the optimal dataset-specific hyperparameters for HDBSCAN clustering.

1. `./03_quast_sourmash_tsne.sh`: it first selects high quality assemblies based on QUAST output (N50>30000, n_contigs<1500, and Total_length between 9000000 & 15000000). Then, creates a sketch of the genomes and calculates mash distances using `sourmash`. The distance matris is used for clustering and subsampling with `t-SNE` and `HDBSCAN`. Finally, it creates samplesheets to perform sensitivity analyses during the schema definition step.

## Step 4 - Gene by gene schema definition

1. `./04_classic_schema_analysis.sh`: as a baseline comparison, we annotate two classic MLST schemas for this parasite (six and three loci respectively) and create phylogenetic trees out of them using multiple sequence alignment as input for `IQTREE`
2. `./05_nfchewbbaca_full.sh`: runs the chewBBACA pipeline with cross-validation for all assemblies that passed QUAST quality control filters
3. `./06_nfchewbbaca_clustered.sh`: runs the sensitivity analyses of the chewBBACA pipeline, using subsets that aim to minimize the weight of isolates from BC

## Step 5 - Downstream analysis

1. `07_interpro_annotation.sh`: conducts annotation of all loci identified in the chewBBACA pipeline for the HQ assemblies as a sanity-check that we are actually identifying previously described giardia proteins

## Helper scripts in bin folder

1. `mlst2dist.py`: chewBBACA companion script to produce table of allelic distances from the tool output.
2. `get_metadata_giardia.py`: companion snippets to obtain sample metadata for all accessions in data set.
2. `quast_analysis.R`: R script to import, clean, and produce summary of assembly QC.
