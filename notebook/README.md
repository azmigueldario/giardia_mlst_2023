# Pipeline analysis v 1.0

## Genome assembly and quality selection

All illumina genomes are assembled using **Shovill** and **Spades**, with default configuration. Using the results from Quality Control (QC) of the assemblies, the contigs below a threshold and overall poor quality draft genomes will be removed before cgMLST analysis.

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
