# Project: **Giardia _cgMLST_ 2023**

## Approach

- Collect publicly available and in-house genomes of Assemblages A an B.
- Normalize and apply quality control tools to raw reads
- Produce annotated and verified genomes for all parasites
- Apply Chewbacca tool to define a core-genome MLST approach

## Datasets

Secondary data analysis project to produce a core-genome multi-locus sequence typing (MLST) approach to _Giardia duodenalis_ assemblages A and B parasites

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

## Usage (v0.1)

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



## Cross validation of cgMLST calling using chewBBACA
