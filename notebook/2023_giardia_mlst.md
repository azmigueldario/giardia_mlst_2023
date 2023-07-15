[toc]

# General information

## Giardia genomes background

- Highly conserved, seems like they are in intense purifying selection.
  Assemblages A and B are the major cause of disease in humans. Others are
  mostly transmitted among animals.  
        + Yet, nucleotide sequence identity between assemblages A and B is around 70%  
        + Assemblages AI and AII have ~99% similarity
- Few introns, mostly localized in genes of tRNA
- HIghly variable genes include:  
		+ Variant-specific Surface proteins (VSP). Found commonly in genomic
		  regions with with rearrangements.  
		+ NEK kinases  
		+ High Cysteine Membrane proteins (HCMP)  
		+ Protein 21.1

## Glossary

- **dN/dS** - Nonsynonymous to synonymous, high when active selection and low
  when purifying selection
- **P15, WB and GS** are commonly used reference strains/genomes. The GS isolate
  is an assemblage B parasite, while the WB strain is the model assemblage A
  organism. 
- **Pseudogenes** are regions with homology to a known gene or a possible ORF
  that lack a functional sequence necessary for transcription or have truncation
  or premature stop codons
- **chewBACCA** is a freely available tool to call gene-by-gene typing schemes.
  Further information can be found in
  https://chewbbaca.readthedocs.io/en/latest/index.html
- **prodigal** is a machine learning algorithm that uses a reference genome for
  an organism as training data to learn how to identify coding sequences (CDS)
  in other query genomes


## Accession numbers for available genomes (Checked on March 16, 2023)

We are looking for short or long read whole genome data of G. duodenalis
assemblages A or B. Initial possible projects include:

1. AACB00000000 / GCA_000002435.2 (Whole  genome of WB isolate (AI))  
2. PRJNA280606 (BCCDC primary study)
3. JXTI00000000 (primary assembly SE reads using 454 technology)
4. [SAMN12878171](https://www.ebi.ac.uk/ena/browser/view/SAMN12878171)  

# Analysis pipeline

## Eagle Environment setup

```sh
# path to scripts for project
SCRIPTS_DIR="/project/60005/mdprieto/giardia_mlst_2023/scripts"

# path to downloaded giardia genomes
BIOR_GENOMES="/project/60005/mdprieto/raw_data/giardia/biorepositories"

# path to BCCDC genomes
BC_GENOMES="/project/60005/mdprieto/raw_data/giardia/BCCDC"

# singularity images
PRODIGAL_IMG="/project/cidgoh-object-storage/images/prodigal_2.6.3.sif"
CHEWBACCA_IMG="/mnt/cidgoh-object-storage/images/chewbacca_3.1.2.sif"
```


## Cedar Environment setup
```sh
```

# Notebook notes

## 20230227 - Preparation

- Created **Pubmed** search to track relevant papers for literature review and
  to fill-up introduction  
        - Finished literature review about genomic studies in Giardia
- Defined a preliminary workflow for QC-Assembly-Annotation in a .ppt file  
        - QC of assembly by backtracking mapping to reference genomes  

## 20230316 - Obtain available relevant genomes

- Searching genomes for _G. duodenalis_ in the European Nucleotide Archive
  (ENA), looking through projects
- SRA search based on organism: _Giardia spp._ and excluding several sequencing
  platforms results in 215 hits. Search details below. 
    + **("giardia"[Organism] AND ("genomic"[Source] OR "other"[Source] OR
      "transcriptomic"[Source])) NOT ("abi solid"[Platform] OR
      "capillary"[Platform] OR "helicos"[Platform] OR "ls454"[Platform]) NOT
      ("amplicon"[Strategy] OR "rna seq"[Strategy])** 
- Downloaded SRR accession list with all available genomes to the project repo,
  contains data from paired end reads of assemblages A and B

## 20230405 - Preparation for chewBACCA run

- With the pipeline for NGS data download of nf-core
  [(nf-core/fetchngs)](https://nf-co.re/fetchngs), I retrieve the genomes from
  the ENA into eagle
- Created a conda environment for the project in Eagle, just in case it is
  necessary
- Then, having the fastq files for all available genomes I separate the metadata
  and raw files into those coming from the local study in the BCCDC and all
  others
- Installed chewbacca singularity image from [depot.galaxy
  link](https://depot.galaxyproject.org/singularity/chewbbaca%3A3.1.2--pyhdfd78af_0)

```sh
# search patterns of a file B in file B, output the not matching lines
grep -v -f fileB.txt fileA.txt > outputFile.txts

# list of BCCDC accessions
BCCDC_ACC="/project/60006/mdprieto/giardia_mlst_2023/processed_data/ACC_BCCDC.txt"

# create a text file with BCCDC data and move it 
for file in $(cat to_move.txt); do mv "$file" /project/60005/mdprieto/raw_data/giardia/BCCDC/fastq; done

#  install chewbacca image
module load apptainer
cd ~/object_images
singularity pull https://depot.galaxyproject.org/singularity/chewbbaca%3A3.1.2--pyhdfd78af_0

```

## 20230412 - Preparation for chewBACCA run

chewBACCA requires training files for prodigal to identify CDS.

- Training files can be produced with prodigal using a reference genome. Thus, I
  use a chromosome-level assembled genome for the WB strain of _G. duodenalis_
  as reference
- Read the complete chewBACCA documentation to get familiar with the tool and
  the necessary dependencies and steps
- Ran preliminary createSchema process

```sh
# downloaded assemblage A reference genome
cd TARGET_DIR
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.gff.gz
```

## 20230515 - Advances

- Performed AlleleCall
- Downloaded **checkm2** singularity image and prepared necessary database
- Created **assembly_qc** script to analyze BCCDC genomes

## 20230620

- After training for a while in nextflow to optimize pipeline, ready to give it
  a go once again
- In a previous meeting, Jimmy Liu suggested to divide randomly
  (cross-validation) my genomes into training and testing 
    - I will create a Nextflow pipeline to quickly do it.
- Also, as fastqc, trimming, trimmed_fastqc, shovill assembly, and assembly
  qc(checkm and quast) is something I do repeatedly, I will implement it in a
  nextflow pipeline
    - I created a draft of the complete pipeline and started developing (check
      CIDGOH QC and nf-seqqc pipeline as template) 
- First step is to unify assembly and QC pipeline for all isolates 
    - Raw data in `fastq` format is already available in my environment as I
      obtained it from NCBI

## 20230621 

- Created an input samplesheet for the pipeline using the command below 

```sh
    # add header line
echo "sample,read1,read2" > input_samplesheet.csv

for read1 in $(ls /home/mdprieto/mdprieto_projects/raw_data/giardia/{BCCDC,repositories}/fastq/*_1.fastq.gz);
    do 
        # get the basename of file and remove suffix
    sample_id=$(echo $read1 | xargs -n 1 basename -s '_1.fastq.gz')
        # replace string '_1' for '_2'
    read2="${read1/_1/_2}"
        # write in a new line for each sample
    echo $sample_id,$read1,$read2 >> input_samplesheet.csv
    done
```

- Started designing nextflow pipeline to run all chewbacca automatically
- **Important note: The selection of training and testing datasets must be done
  outside nextflow as it is a random process**. Otherwise, nextflow will not be
  able to resume. 
- Division intro training and testing datasets will be done with Python, using
  the scikit-learn module. We can prepare an environment for this process as
  follows:
  
```sh
conda create -n sklearn_env -c conda-forge scikit-learn pandas
conda activate sklearn_env
```

- Created first instance of python script to divide dataset into test/training,
  and do it 'n' times according to cross-validation
    - Got more familiar with the usage of `argparse` module to parse arguments
      provided to the script in the command line
- Drafted preliminary `main.nf` for chewBACCA pipeline

## 20230711 - Creating input channels for cross-validation

- Different ideas to classify reads as input channels
    - Have a long format with columns {...set1,set2,set3,...setn} and divide
      them by hardcoding what column will be kept as train/test identifier
    - Use a long format and just filter based on the value contained in the set
      column
    - Create a custom groovy script that maps the samples however I want them

## 20230713 - Developing workflow

- Created new samplesheet that contains contigs only, that is necessary as input for chewbacca

```sh
echo "sample,contig" > input_samplesheet_contig.csv

for contig in $(ls /home/mdprieto/mdprieto_projects/raw_data/giardia/{BCCDC,repositories}/fasta/*fasta);
    do 
        # get the basename of file and remove suffix
    sample_id=$(echo $contig | xargs -n 1 basename -s '_FINAL.contigs.polished.fasta')
        # write in a new line for each sample
    echo $sample_id,$contig >> input_samplesheet_contig.csv
    done
```
