#!/bin/bash                                 
#SBATCH --mem-per-cpu=6G                    
#SBATCH --time=09:30:00                     
#SBATCH --cpus-per-task=6                  
#SBATCH --job-name="sra_download_giardia"           
#SBATCH --chdir=/scratch/mdprieto/          
#SBATCH --output=jobs_output/giardia_assembly/%x_%j.out 

#########################################################################################################
#                        Preparation 
#########################################################################################################

# load necessary modules
module load apptainer nextflow

# define INPUT list and OUTPUT directory
project_repo="/project/60006/mdprieto/giardia_mlst_2023/"
accessions="${project_repo}/processed_data/accessions"
insdc_genomes="${project_repo}/input_data/giardia/test"
outdir_ref="${insdc_genomes}/../assemblage_A"

# path to sra-tools image
sra_tools_container="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-sra-tools-3.2.1--h4304569_1.sif"

# merge all accessions into single file
cat ${accessions}/BCCDC_accessions.txt ${accessions}/SRR_accessions_2025.csv > ${accessions}/all_accessions_2025.csv

# ncbi datasets image
sra_tools_img="project/60006/mdprieto/giardia_mlst_2023/sra-tools_latest.sif "

# make directories if necessary
mkdir -p ${insdc_genomes} ${outdir_ref}

#########################################################################################################
#                           Download paired-end reads using SRA tools
#########################################################################################################

# download fastq data
while read -r SRR; do
    echo "Downloading $SRR..."
    
    apptainer exec ${sra_tools_container} \
        fasterq-dump \
        --split-files \
        --outdir raw_fastqs \
        --threads 6 \
        --progress \
        $SRR
        
done < ${accessions}/all_accessions_2025.csv

#########################################################################################################
#                               Reference genome
#########################################################################################################

# download chromosome level assembly of giardia Duodenalis (Illumina + PacBio)
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz \
    --output  ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.fna
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.gff.gz \
    --output  ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.gff

# leave a copy of ref genome (for pipeline) in a fasta subdirectory of the repository accessions
mkdir -p ${insdc_genomes}/wb_reference && \
    cp ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.fna $(insdc_genomes)/wb_reference