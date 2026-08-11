#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description: Download reference fasta(s) for specific Giardia loci included in 
#              previously published MLST schemas
# Last modified: 2026-03-10
# Usage: ./quast_sourmash_tsne.sh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load configuration file and create output folder(s)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

source "./../../config_paths.sh"
mkdir "${input_dir}/classic_schemas"
mkdir "${input_dir}/classic_schemas/six_loci"
mkdir "${input_dir}/classic_schemas/three_loci"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Download reference fasta(s) for each loci in the schemas
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

curl -X GET "https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/5698276,5699725,5700156,5697398,5699357,5697828,5701000,5699927,5700811/download?include_annotation_type=FASTA_PROTEIN&include_annotation_type=FASTA_GENE" \
 -H 'accept: application/zip' \
 --output "${input_dir}/classic_schemas/MLST.zip"

# Unpack only fasta file
unzip \
    -j "${input_dir}/classic_schemas/MLST.zip" "*.fna" \
    -d "${input_dir}/classic_schemas"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load seqkit and separate multifasta file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

module load StdEnv/2023 seqkit/2.5.1

# Clean names and split fasta
seqkit replace \
    --pattern "\.1:[c0-9]+-[0-9]+ " \
    --replacement "-1:" \
    "${input_dir}/classic_schemas/gene.fna" |
seqkit split \
    --by-id \
    --line-width 80  \
    --id-regexp "(^.*GL[0-9_]+) " \
    --two-pass \
    --out-dir "${input_dir}/classic_schemas"

cd "${input_dir}/classic_schemas" &&
rename "stdin.part_" "" *

# Separate by schema
mv *GL50803_004812* *GL50803_0093938* *GL50803_0021942* \
    "${input_dir}/classic_schemas/three_loci"

mv *.fasta "${input_dir}/classic_schemas/six_loci"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Download outgroup assembly (Giardia muris)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

wget \
    -O  "${input_dir}/classic_schemas/GCA_006247105-1_1_genomic.fna.gz" \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/006/247/105/GCA_006247105.1_UU_GM_1.1/GCA_006247105.1_UU_GM_1.1_genomic.fna.gz

# Extract fasta
gunzip \
    -d "${input_dir}/classic_schemas/GCA_006247105-1_1_genomic.fna.gz" > "${input_dir}/classic_schemas/GCA_006247105-1_1_genomic.fna"
