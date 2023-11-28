# Giardia MLST creation and validation

## General information

### Background

- Highly conserved genome, seems like they are in intense purifying selection.
  Assemblages A and B are the major cause of disease in humans. Others are
  mostly transmitted among animals.  
  - Yet, nucleotide sequence identity between assemblages A and B is around 70%.
  Assemblages AI and AII have ~99% similarity.
- Few introns, mostly localized in genes of tRNA, so it behaves somehow like a bacterial genome.
- HIghly variable genes include:  
  - Variant-specific Surface proteins (VSP). Found commonly in genomic
    regions with with rearrangements.  
  - NEK kinases  
  - High Cysteine Membrane proteins (HCMP)  
  - Protein 21.1

### Objectives

1. Find and retrieve all available genomes for Giardia intestinalis (Assemblages A and B)
2. Create a robust core-genome MLST scheme for future epidemiological analysis of outbreaks
    - Cross-validation of results given the small number of available datasets

### Repository/folder structure

**_main_** contains tested scripts and results, will contain final pipeline
**_dev_** is mainly for developing purposes and testing in eagle
**_cedar_** is a local branch used in Compute Canada (Digital Research Alliance) for analysis

To not overwrite local changes in **_cedar_** I will use the following snippet to update it:

```sh
git stash
git merge <branch>
git 
```

As usual, the project is divided in four main subdirectories and inside a github repository (<https://github.com/azmigueldario/giardia_mlst_2023>)

- **notebook:** contains the _markdown_ file that documents all the advances and troubleshooting. It also contains a graph summarizing the analytical pipeline
- **output:** contains all results from analysis
- **processed_data:** is where all input samplesheets, input datasets, and relevant information to run the pipeline are located.
  - **NOTE:** The _raw data_ for this project is saved outside the repository.
- **scripts:** contains all analytical scripts and workflows developed for this project

```sh
.
├── notebook
├── output
├── processed_data
│   └── cross_validation_input
└── scripts
    └── nextflow
```

### Glossary

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
  <https://chewbbaca.readthedocs.io/en/latest/index.html>
- **prodigal** is a machine learning algorithm that uses a reference genome for
  an organism as training data to learn how to identify coding sequences (CDS)
  in other query genomes
- **cross-validation** using the same set of data to evaluate accuracy of a model by partitioning into training and testing datasets and conducting analysis iteratively

### Accessions of available data (verified 20230316)

We are looking for short or long read whole genome data of G. duodenalis
assemblages A or B. Initial possible projects include:

1. AACB00000000 / GCA_000002435.2 (Whole  genome of WB isolate (AI))  - **Reference genome**
2. PRJNA280606 (BCCDC primary study)
3. JXTI00000000 (primary assembly SE reads using 454 technology)
4. [SAMN12878171](https://www.ebi.ac.uk/ena/browser/view/SAMN12878171)  

### Environment setup

```sh
#################################### Eagle ####################################

# path to scripts for project
SCRIPTS_DIR="/project/60005/mdprieto/giardia_mlst_2023/scripts"

# request interactive session for development
salloc --time=02:00:00 --ntasks=1 --cpus-per-task=8  --mem-per-cpu=10G 

#################################### Cedar ####################################

# path to giardia scripts directory
SCRIPTS_GIARDIA="/project/6056895/mdprieto/giardia_mlst_2023/scripts"

# interactive allocation
salloc --time=2:00:0 --ntasks=1 --cpus-per-task=8  --mem-per-cpu=8G --account=def-whsiao-ab
```

## Notebook notes

### 20230227 - Preparation

- Created **Pubmed** search to track relevant papers for literature review and
  to fill-up introduction  
        - Finished literature review about genomic studies in Giardia
- Defined a preliminary workflow for QC-Assembly-Annotation in a .ppt file  
        - QC of assembly by backtracking mapping to reference genomes  

### 20230316 - Obtain available relevant genomes

- Searching genomes for _G. duodenalis_ in the European Nucleotide Archive
  (ENA), looking through projects
- SRA search based on organism: _Giardia spp._ and excluding several sequencing
  platforms results in 215 hits. Search details below.
  - **("giardia"[Organism] AND ("genomic"[Source] OR "other"[Source] OR
      "transcriptomic"[Source])) NOT ("abi solid"[Platform] OR
      "capillary"[Platform] OR "helicos"[Platform] OR "ls454"[Platform]) NOT
      ("amplicon"[Strategy] OR "rna seq"[Strategy])**
- Downloaded SRR accession list with all available genomes to the project repo,
  contains data from paired end reads of assemblages A and B

### 20230405 - Preparation for chewBACCA run

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

### 20230412 - Preparation for chewBACCA run

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

### 20230515 - Advances

- Performed AlleleCall
- Downloaded **checkm2** singularity image and prepared necessary database
- Created **assembly_qc** script to analyze BCCDC genomes

### 20230620

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

### 20230621

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

### 20230711 - Creating input channels for cross-validation

- Different ideas to classify reads as input channels
  - Have a long format with columns {...set1,set2,set3,...setn} and divide
      them by hardcoding what column will be kept as train/test identifier
  - Use a long format and just filter based on the value contained in the set
      column
  - Create a custom groovy script that maps the samples however I want them

### 20230713 - Developing workflow

- Created new samplesheet that contains contigs only, that is necessary as input for chewbacca
- Prepared question about how to approach cross-validation for nextflow forum. I have thought about several ideas but have no final solution and do not want to stale advances anymore.
  - Includes a sample_sheet in long or wide format, any approach would work and a sample of processes to run.
- Created a draft of processes to run in ChewBACCA

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

### 20230717 - Developing workflow: adjusting python script for test/sample split

Modified python script to produce consolidated long and wide format `.csv`. The new script produces results in long or wide format too.

- Used `pandas` library to handle data frames in comma separated files
  - `pd.assign` lets me add additional columns of variables
  - `dataset.loc[dataset[var]==condition]` is helpful to evaluate a condition and modify a value inside the dataframe
  - For string interpolation in **Python** I have to remember to use `f'string{replacement}'`
- A new pilot dataset is created using head of the main dataset `df.head(20)` and then doing the random subsampling
- Still working on setting up the right way to feed all data for cross-validation into ChewBACCA
- Created a `input_channels.nf` file to test the pipeline in a subset of the data

```nf
// Using branch and setting according to match in set column

Channel
    .fromPath(params.input_file)
    .splitCsv(header: true)
    .branch { row -> 
            set1: row.set == "set1"
                tuple(row.sample, row.set, row.value, row.contig) 
            set2: row.set == "set2"
                tuple(row.sample, row.set, row.value, row.contig) 
            set3: row.set == "set3"
                tuple(row.sample, row.set, row.value, row.contig) 
                }
    .set{fasta_ch}

```

#### Split csv file in number of lines but keeping header

```sh
tail -n +2 INPUT_FILE | split -l 4 - prefix_
for file in prefix_*
do
    head -n 1 INPUT_FILE > tmp_file
    cat $file >> tmp_file
    mv -f tmp_file $file
done
```

### 20230717 - Developing workflow: setup initial process for prodigal training

- When setting up a pipeline with **Singularity** in nextflow, I have to add the following parameters to the `nextflow.config` file to create the mounting environment to execute the code and to load singularity
    > singularity.enabled = true
    > singularity.autoMounts = true
- Defined process to create prodigal training file
- Decided to create a working pipeline with a single csv input and repeat it ten times. Once I am more proficient with nextflow I can come back and improve upon it.
- For the CreateSchema and probably other processes, I need to `collect().flatten()` the contigs to feed them and I still do not know how to do it
  - Singularity image is not pulling correctly from '<https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0>'. The following code was successful in a pilot run:

```sh
singularity pull --name depot.galaxyproject.org-singularity-chewbbaca%3A3.2.0--pyhdfd78af_0.img https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0\
```

### 20230726

- Added output directory to process, optimized writing of first two processes
- Improved nextflow.config to parametrize resources for each process
- Updated definitions of output
- Currently testing AlleleCall and RemoveParalogs processes

### 20230731

- Tried to allow input for several sets inside the same `.csv` samplesheet using the `multimap()` groovy operator. However,
- ALLELE_CALL process is working adequately now, may require matching to select final output to copy in results folder
- Trying to use bactopia to produce standardized assemblies of available data
  - Working pretty slowly as nextflow distributes jobs only in the CPUS of a single node and 8 CPUS are available per node

```sh
# create bactopia tab separated samplesheet
echo -e "sample\truntype\tr1\tr2\textra" > bactopia_samplesheet.csv

for read1 in $(ls /project/6056895/mdprieto/raw_data/giardia/{BCCDC,repositories}/fastq/*_1.fastq.gz);
    do
        # get the basename of file and remove suffix
    sample_id=$(echo $read1 | xargs -n 1 basename -s '_1.fastq.gz')
        # replace string '_1' for '_2'
    read2="${read1/_1/_2}"
        # write in a new line for each sample
    echo -e "$sample_id\tpaired-end\t$read1\t$read2"  >> bactopia_samplesheet.csv
    done 
```

### 20230802 - Developing chewBBACA pipeline, troubleshoot bactopia in cedar

- Updated `draft.nf` now it includes a working copy of `remove_paralogs` and `extract_cgmlst` processes
- Improved naming of input/output throughout the pipeline to better reflect the usage of train/test subsets
- Troubleshooting execution of bactopia in **cedar** cluster, import error for `MINMER_SKETCH` process
  
  - Seems to be an issue where the container is using the `python` interpreter and `numpy` module from the cedar environment, could be solved by adding an option to isolate the container environment to the bactopia `docker.config` (inherited to `singularity.config`)

```nf
  # worked
withLabel: 'minmer_sketch|minmer_query' {
    container = "${docker_repo}/minmers:${manifest.version}"
    singularity.runOptions = "-c --bind /project,/scratch"
    }
```

### 20230804 - Assembling all Illumina genomes in batch (Cedar)

- Now, there is an issue where `spades` is running out of memory for its tasks. I modified the `base.config` file of Bactopia to add more memory for `assemble_genome` process
  - Additional memory did not work, seems like the HPC `/tmp` folder may not have sufficient memory for `Spades` so I try adding an additional option for the **Shovill** assembler in bactopia: `--shovill_opts "--tmp-dir /scratch/mdprieto/tmp"`
  - To speed up the assembly and other processes, I also increase the maximum runtime and the cpus of certain tasks in bactopia. It is necessary as the Giardia genome may be significantly larger than a bacterial one (12M vs 6.5M)

### 20230805 - Assembling all Illumina genomes in batch (Cedar, Eagle)

- Taking too long to schedule a long running job in Cedar
- Created local config file for execution of **Bactopia** through slurm submission of jobs in Eagle (`scripts/eagle.config`) that includes:
  - runOptions to map a large directory to the `/tmp` inside singularity containers
  - detailed instructions to run jobs inside slurm on Eagle
  - configuration for processes to run for a longer time and not be cancelled due to large size of Giardia spp. genome (compared to bacteria)

### 20230830 - restart work on the project

- Rebased the development branch with the main branch for this project
- Input datasets were removed by mistake, so I recreate them. Yet, there is poor documentation in the python script so I improve it for future uses including sanity check for input dataset as well as detailed explanations of the expected parameters.
  - The ideal output is a long_format file that has four columns containing [sample, contig, set, value] where the value holds the assignment of train or test in every iteration.

### 20230831 - Continue development of pipeline

- In order to process several sample sets using the same prodigal training file, the latter must be assigned to a value channel. Thus, I use the `first()` operator as follows:
  > prodigal_ch = PRODIGAL_TRAINING().first()
  
### 20231012 - Preparing all assembly inputs

NCBI accessions:

- I cleaned the list of NCBI accessions to remove the reference genome for Giardia Assemblage A. This one is downloaded directly using `curl`. It is now ready for prediction with prodigal and also to be included in the MLST pipeline
- Updated all pathways in eagle as the folder structure has changed in the last couple of months. Also, no need to assemble reference genome, as it can be downloaded after polishing.
- Updated version (1.10.1) of `fetchngs` is now working, described all steps to reproduce analysis in `./notebook/README.md`
  - Does not require `--force_sra_tools` option

```sh
# verifies data download script
sbatch /project/60006/mdprieto/giardia_mlst_2023/scripts/download_data_repositories.sh
```

### 20231116 - Catch up with progress

- Read all previous documentation, commited pending changes, and verified download of all datasets
- I am verifying the assembly of all genomes (short and long reads) and the long-reads assembler is running out of memory (will increase to at least 32)
  - Bactopia adjusts maximum memory to the one available in the primary job (the one that submits the other processes), so I set up `-profile slurm` to avoid that
  - Also increased the maximum amount of memory to 48 for resubmitted jobs and standard runtime to 6h for assembly (was hanging)
    - Must specify the target queue in eagle to avoid errors with the default parameters `--slurm_queue 'cpubase_bycore_b1'`
      - NOTE: Above is not necessary as long as the `queue = 'cpubase_bycore_b1'` is added inside the process configuration for the slurm profile

### 20231127 - Troubleshoot assembly pipeline

- Implemented a new profile in bactopia.config that executes in slurm and increases resources for assembly while limiting the number of jobs submitted every minute.
  - Pilot Bactopia running **WORKING**
- Added bash script to produce bactopia samplesheet (`scripts/bin`) when needed

### 20240116 - Reran pipeline

- By mistake, the original results of assembly of all data were deleted. As no backup is available, I run the job again (`/scripts/assembly_and_qc.sh`)

## Progress and pending work

### Completed

- Downloaded all Giardia Assemblage A and B genomes (Verified)
- Produced standardized fasta assemblies using SPAdes (Verified)

### Testing example

```sh
SAMPLESHEET="/project/60006/mdprieto/giardia_mlst_2023/processed_data/pilot_bactopia_samplesheeet.csv"
CUSTOM_CONFIG="/project/60006/mdprieto/giardia_mlst_2023/scripts/eagle_bactopia.config "

###############################################################################################

nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    -resume \
    --nfconfig $CUSTOM_CONFIG \
    --samples $SAMPLESHEET \
    --outdir /scratch/mdprieto/results/bactopia_giardia \
    --shovill_assembler spades \
    --skip_amr \
    --long_reads \
    --skip-prokka \
    --skip_mlst \
    --max_memory 90 \
    --cleanup_workdir
```
