# Project: **Giardia _cgMLST_ 2023**

## Approach

- Collected publicly available and in-house sequenced genomes of Assemblages A an B (`sra-tools`)
- Genome assembly and quality control using the [Bactopia](https://bactopia.github.io/v3.0.0/) pipeline v3.0.0
    - Normalized and performed quality control in raw reads to prune low quality input  (`fastQC + multiQC`)
    - Independent runs for Illumina short reads and hybrid (Illumina + ONT) assemblies  
- Define assemblage of parasites by competitive mapping of draft assemblies against references of both assemblages (`samtools`)
- Create nextflow pipeline with `chewBBACA` tool commands to define a core-genome MLST schema and evalute it in a minimum spanning tree (`grapetree`)
- Sensitivity analyses: 
    - subsampled input genomes to reduce overrepresentation of BC, Canada genomes and re-ran nextflow chewBBACA pipeline (`mash + t-SNE + HDBSCAN`)
    - created phylogenetic trees for the same input genomes used in the cgMLST schema using classic 3 or 6-loci schemas (`chewBBACA + MAFFT + IQTREE`)

## Datasets

Secondary data analysis of Giardia intestinalis assembles A and B (parasites of human interest) to produce a core-genome multi-locus sequence typing (MLST) schema. The sequencing files were obtained through a comprehensive search of the INSDC in November 2024 (`PRJNA561185`, `PRJNA280606` and `PRJEB3213`) and an update in December 2025 (`PRJNA1110996`)

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

## Usage (v0.4)

### Preparation

### 1. Data retrieval and cleaning of metadata

Download data and sample-level metadata using the `sra_tools` software. A custom output path can be specified by modifying the variable `insdc_genomes`.
A simple `curl` command was used to obtain reference sequences and reference annotation from the NCBI accession **GCF_000002435.2** and from another set of high quality hybrid assemblies combining PacBio + Illumina

```sh
./scripts/data_retrieval_assembly/download_data_repositories.sh
```

## 2. Genome assembly and quality selection

For this project, we run **Bactopia v3.0.0** and use Apptainer/Singularity containers to run the commands in a HPC. 

Before creating the input samplesheet, we sanitized the filenames to remove low hypens "_" as they produce conflicts with the parsing logic of **Bactopia**.

- To create the samplesheet it is recommended to run the command below, which is available in the conda environment of the **Bactopia** tool. If **Conda** is not available, you can use a Bactopia container image to run the command.

```sh
bactopia_img="/path/to/bactopia/apptainer.sif"
apptainer exec --cleanenv ${bactopia_img} \
    bactopia prepare \
        --path $insdc_genomes/fastq/
```

All illumina genomes are assembled using **Shovill** and **Spades**, with default configuration, inside the **Bactopia** nextflow pipeline. 

- Using the results from the assembly Quality Control (QC) using **QUAST**, the contigs below a threshold and overall poor quality draft genomes will be removed before cgMLST analysis. Samples with hybrid assembly are independently processed and assembled with **Bactopia**. 

```sh
# Illumina only reads
./scripts/data_retrieval_assembly/assembly_and_qc.sh

# Samples with hybrid sequencing
./scripts/data_retrieval_assembly/hybrid_assembly.sh
```

## 3. Filter datasets cross-validation

First, we review the results of the assembly process and prune the samples that have low N50 (<30,000), a large number of contigs (n > 1300), or a genome size outside the range of the reference assembly (+/- 20%)

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
