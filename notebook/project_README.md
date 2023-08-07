# Pipeline analysis v 1.0

## Genome assembly and quality selection

All illumina genomes are assembled using **Shovill** and **Spades**, with default configuration. Using the results from Quality Control (QC) of the assemblies, the contigs below a threshold and overall poor quality draft genomes will be removed before cgMLST analysis.

```sh
# readsQC, assembly, assemblyQC
scripts/cedar/full_assembly_and_qc.sh
```

## Cross validation of cgMLST calling using chewBBACA
