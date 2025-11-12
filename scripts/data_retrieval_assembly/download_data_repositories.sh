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
project_root="/project/60006/mdprieto/giardia_mlst_2023/"
accessions="${project_root}/processed_data/accessions"
insdc_genomes="${project_root}/input_data/giardia/repositories"
outdir_ref="/mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A"

# merge all accessions into single file
cat $accessions/ACC_BCCDC.txt $accessions/SRR_ACC_list.txt > $accessions/ALL_ACC.csv

####################################  Download #####################################################

# download fastq data
nextflow run nf-core/fetchngs -r 1.10.1 \
    --input "${accessions}/ALL_ACC.csv" \
    --outdir ${insdc_genomes} \
    -profile singularity \
    --nf_core_pipeline 'taxprofiler' \
    -resume

################################  Reference genome ###############################################

# download chromosome level assembly of giardia Duodenalis (Illumina + PacBio)
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz \
    --output  ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.fna
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.gff.gz \
    --output  ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.gff

# leave a copy of ref genome (for pipeline) in a fasta subdirectory of the repository accessions
mkdir -p ${insdc_genomes}/wb_reference && \
    cp ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.fna $(insdc_genomes)/wb_reference