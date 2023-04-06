[toc]

# General information

## Giardia genomes

- Highly conserved, seems like they are in intense purifying selection. Assemblages A and B are the major cause of disease in humans. Others are mostly transmitted among animals.  
        + Yet, nucleotide sequence identity between assemblages A and B is around 70%  
        + Assemblages AI and AII have ~99% similarity
- Few introns, mostly localized in genes of tRNA
- HIghly variable genes include:  
		+ Variant-specific Surface proteins (VSP). Found commonly in genomic regions with with rearrangements.  
		+ NEK kinases  
		+ High Cysteine Membrane proteins (HCMP)  
		+ Protein 21.1

## Glossary

- **dN/dS** - Nonsynonymous to synonymous, high when active selection and low when purifying selection
- **P15, WB and GS** are commonly used reference strains/genomes. The GS isolate is an assemblage B parasite. 
- **Pseudogenes** are regions with homology to a known gene or a possible ORF that lack a functional sequence necessary for transcription or have truncation or premature stop codons


## Accession numbers for available genomes (Checked on March 16, 2023)

We are looking for short or long read whole genome data of G. duodenalis assemblages A or B. 

1. AACB00000000 / GCA_000002435.2 (Whole  genome of WB isolate (AI))  
2. PRJNA280606 (BCCDC primary study)
3. JXTI00000000 (primary assembly SE reads using 454 technology)
4. [SAMN12878171](https://www.ebi.ac.uk/ena/browser/view/SAMN12878171)  

# Analysis pipeline

## Environment setup - Eagle
```sh
# path to scripts for project
SCRIPTS_DIR="/project/60005/mdprieto/giardia_mlst_2023/scripts"
```

## Environment setup - Cedar
```sh
```

# 20230227

- Created **Pubmed** search to track relevant papers for literature review and to fill-up introduction  
        - Finished literature review about genomic studies in Giardia
- Defined a preliminary workflow for QC-Assembly-Annotation in a .ppt file  
        - QC of assembly by backtracking mapping to reference genomes  
# 20230316 

- Searching genomes for _G. duodenalis_ in the European Nucleotide Archive, looking through projects
- SRA search based on organism: _Giardia spp._ and excluding several sequencing platforms results in 215 hits. Search details below. 
> ("giardia"[Organism] AND ("genomic"[Source] OR "other"[Source] OR "transcriptomic"[Source])) NOT ("abi solid"[Platform] OR "capillary"[Platform] OR "helicos"[Platform] OR "ls454"[Platform]) NOT ("amplicon"[Strategy] OR "rna seq"[Strategy]) 
- Downloaded SRR accession list with all available genomes to the project repo


