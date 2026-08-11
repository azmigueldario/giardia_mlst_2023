### Genomic subtyping schema for Giardia duodenalis 

## Introduction

1. Global epidemiology of Giardia duodenalis
    - Classification of assemblages
    - Clinical burden for humans and outbreaks
2. Genomics and surveillance
    - Genomic characteristics: ploidy, classification, virulence 
    - Microbiological surveillance: available tools (culture and typification schemas)
    - Opportunity for direct sequencing and need for genomic subtyping (gap)

## Methods

1. Dataset definition (n = 180)
    - Local collection of Giardia genomes (ILLUMINA)
    - Enriched with publicly available genomes
2. Quality control for input data 
    - Competitive mapping against reference assemblages (A and B) affecting humans
    - K-mer classification for contaminated samples
    - Metadata description (Harmonization)
    - Standardized read filtering and assembly with SPADES (Bactopia)
    - Filter out assemblies with poor quality 
3. cgMLST definition
    - chewBBACA nextflow pipeline for cross-validation and cgMLST calling (venn diagram)
    - Sensitivity analyses for cgMLST definition:
        - Subsample of available specimens after clustering (maximize sketching distances)
        - Subsample of available specimens after clustering (random)
    - Annotation of selected loci
    - Minimum spanning trees for each schema dedi
        - Comparison against phylogenetic tree

## Results

1. Summary of dataset: project accessions, sample size, and brief description of dataset (sequencing methods, isolation source, assemblages)
    - Table 1 (Characteristics of data used)
2. Classic MLST schemas: resolution, definition of loci used
    - Figure 2 (ML trees of classic MLST schemas)
3. Novel gene by gene schema: representativeness of data
    - Figure 3:
        - Summary of pipeline for definition of gene by gene schema with cross-validation
        - Minimum spanning tree of Giardia isolates using the novel schema
    - Figure 4: 
        - Sensitivity analyses of gene-by-gene schema reducing the overrepresentation of BC samples (2 MST trees)

## Discussion


## Supplementary info

1. Table 1: annotation of loci in gene by gene schema
 