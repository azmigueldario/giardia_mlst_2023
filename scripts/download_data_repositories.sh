#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                    
#SBATCH --time=09:30:00                     
#SBATCH --cpus-per-task=6                  
#SBATCH --job-name="sra_download_giardia"           
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/%x_%j.out 

################################### Preparation #####################################################

# load necessary modules
module load apptainer nextflow

# define INPUT list and OUTPUT directory
ACCESIONS="/project/60006/mdprieto/giardia_mlst_2023/processed_data/accessions"
BIOR_GENOMES="/project/60006/mdprieto/raw_data/giardia/"

# merge all accessions into single file
cat $ACCESIONS/ACC_BCCDC.txt $ACCESIONS/SRR_ACC_list.txt > $ACCESIONS/ALL_ACC.csv

####################################  Download #####################################################

# download fastq data
nextflow run nf-core/fetchngs -r 1.10.1 \
    --input "$ACCESIONS/ALL_ACC.csv" \
    --outdir $BIOR_GENOMES \
    -profile singularity \
    --nf_core_pipeline 'taxprofiler' \
    -resume

################################  Reference genome ###############################################

# download chromosome level assembly of giardia Duodenalis (Illumina + PacBio)
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz \
    --output  /mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.fna
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.gff.gz \
    --output  /mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.gff

# leave a copy of ref genome (for pipeline) in a fasta subdirectory of the repository accessions
mkdir -p ${BIOR_GENOMES}/fasta && \
    cp /mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.fna $BIOR_GENOMES/fasta/