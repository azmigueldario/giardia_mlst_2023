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
project_repo="/project/60006/mdprieto/giardia_mlst_2023/"
accessions="${project_repo}/processed_data/accessions"
insdc_genomes="${project_repo}/input_data/giardia/repositories"
outdir_ref="${insdc_genomes}/../assemblage_A"

# merge all accessions into single file
cat ${accessions}/ACC_BCCDC.txt ${accessions}/SRR_ACC_list.txt > ${accessions}/ALL_ACC.csv

# ncbi datasets image
ncbi_datasets_img="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-ncbi-datasets-cli:14.26.0.sif"

##########################################  Download fetch NGS - Illumina #####################################

# download fastq data
nextflow run nf-core/fetchngs -r 1.10.1 \
    --input "${accessions}/ALL_ACC.csv" \
    --outdir ${insdc_genomes} \
    -profile singularity \
    --nf_core_pipeline 'taxprofiler' \
    -resume

################################ FTP - Extra PacBio assemblies #######################################

cat <<EOF > /tmp/assemblies.txt
GCA_029168855.1
GCA_029168835.1
GCA_029168845.1
GCA_029168795.1
GCA_029168785.1
GCA_029168805.1
GCA_029168735.1
GCA_029168715.1
GCA_029168725.1
EOF

# Create a folder for the files
mkdir -p fastas
while read acc; do
    prefix=$(echo $acc | sed 's/...$//' | sed 's/.\{3\}/&\//g')
    url="https://ftp.ncbi.nlm.nih.gov/genomes/all/${prefix}${acc}/"
    echo "Downloading $acc from $url"
    wget -r -np -nH --cut-dirs=7 -A "*genomic.fna.gz" "$url"
done < /tmp/assemblies.txt

################################  Reference genome ###############################################

# download chromosome level assembly of giardia Duodenalis (Illumina + PacBio)
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.fna.gz \
    --output  ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.fna
curl https://ftp.ncbi.nlm.nih.gov/genomes/refseq/protozoa/Giardia_intestinalis/latest_assembly_versions/GCF_000002435.2_UU_WB_2.1/GCF_000002435.2_UU_WB_2.1_genomic.gff.gz \
    --output  ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.gff

# leave a copy of ref genome (for pipeline) in a fasta subdirectory of the repository accessions
mkdir -p ${insdc_genomes}/wb_reference && \
    cp ${outdir_ref}/GCF_000002435.2_UU_WB_2.1_genomic.fna $(insdc_genomes)/wb_reference