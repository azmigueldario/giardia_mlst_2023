# Project: **Giardia _cgMLST_ 2023**

## Approach

- Collected publicly available and in-house sequenced genomes of Assemblages A an B (`nfcore-fetchngs`)
- Normalized and performed quality control in raw reads to prune fragmented and low quality assemblies (`fastQC + multiQC`)
- Define assemblage of parasites by comprosite (`samtools`)
- Create nextflow pipeline with `chewBBACA` tool commands to define a core-genome MLST schema and evalute it
- Sensitivity analysis 1: subsampled input genomes to reduce overrepresentation of BC, Canada genomes and re-ran nextflow chewBBACA pipeline (`mash + t-SNE + HDBSCAN`)
- Sensitivity analysis 2: created phylogenetic trees for the same input genomes used in the cgMLST schema using classic 3 or 6-loci schemas (`chewBBACA + MAFFT + IQTREE`)

## Datasets

Secondary data analysis of Giardia intestinalis assembles A and B (parasites of human interest) to produce a core-genome multi-locus sequence typing (MLST) schema

## Analysis steps

1. Obtain the datasets
2. Standardized approach to assembly of the datasets and assembly QC
3. Calculate genomic distances using `mash` sketches
4. Filter out some assemblies based on quality QC and clustering **[t-sne and HDBSCAN]**
5. Run nextflow pipeline **[chewBBACA + crossvalidation]**
6. Obtain hamming distances of cgMLST calls
7. Evaluate annotation and quality of schema
8. Produce minimum spanning tree

## Repository organization

Every analysis folder contains, if applicable, a `README.md` file and instructions to reproduce the computing environment. 

```sh
.
├── notebook
├── output
│   ├── genomes_list
│   ├── nf_chewbbaca
│   ├── quast
│   └── smash
├── processed_data
│   ├── accessions
│   ├── bactopia_samplesheets
│   └── metadata
└── scripts
    ├── classic_schemas
    ├── clustering_tsne
    ├── data_processing_assembly
    ├── nf_chewbbaca
    └── plots_dowstream_analysis
```

## Usage (v0.2)

## Data retrieval and cleaning of metadata

Download data and sample-level metadata using the `nfcore-fetchngs` pipeline.
Curl was used to obtain reference sequences and reference annotation.

```sh
./scripts/data_processing_assembly/download_data_repositories.sh
```



## Genome assembly and quality selection

All illumina genomes are assembled using **Shovill** and **Spades**, with default configuration. Using the results from Quality Control (QC) of the assemblies, the contigs below a threshold and overall poor quality draft genomes will be removed before cgMLST analysis. Samples with hybrid assembly are also processed by bactopia and typically yield better results. 

- Download all raw sequences from NCBI. The script requires nextflow and singularity (Apptainer). It uses NCBI tools, which are notoriously inconsistent, so try a few times if necessary.
- Everything will be downloaded to the PATH you specify as `BIOR_GENOMES`

```sh
./scripts/download_data_repositories.sh
```

- Create a samplesheet for input to the assembly and QC pipeline (Bactopia). We use the python package that complements bactopia, made available through a [Singularity image of Bactopia v3.0](https://depot.galaxyproject.org/singularity/bactopia%3A3.0.0--hdfd78af_0).
    1. Some files may have name conventions that make it harder for the script to read, so make sure that there are no low hyphens "_"

```sh
# readsQC, assembly, assemblyQC
scripts/cedar/full_assembly_and_qc.sh
```

## Prepare datasets for cross-validation

### Clustering and subsampling

## Cross validation of cgMLST calling using nf-chewBBACA

### Calculate final cgMLST distances

```sh
# HQ samples
python3 scripts/python/6_mlst2dist.py \
    processed_data/nf_chewbbaca/full/cgMLST90.tsv \
    processed_data/nf_chewbbaca/full/cgMLST90_distances.tsv \
    --outfmt 'TSV'

# clustered samples
python3 scripts/python/6_mlst2dist.py \
    processed_data/nf_chewbbaca/clustered/cgMLST90.tsv \
    processed_data/nf_chewbbaca/clustered/cgMLST90_distances.tsv \
    --outfmt 'TSV'
```

### Minimum spanning tree

The `grapetree` tool is ideal to produce and manipulate the plot parameters of a minimum spanning tree.

- To reproduce the conda environment for this tool, use the requirements file: `scripts/python/grapetree.yml`
- The IDE is started by running the `grapetree` command and takes as input:
  - A profile tab-delimited file (like `processed_data/nf_chewbbaca/clustered/cgMLST90.tsv`)
  - A metadata file to customize the MST presentation (`processed_data/metadata/giardia_metadata.tsv`)

Our plots have branches in Log-scale to avoid overcrowding, are colored by sampling location, and collapse into a single node samples with less than 60 alleles of separation.
