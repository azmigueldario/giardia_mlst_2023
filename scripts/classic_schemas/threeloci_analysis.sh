#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH --time=00:04:00
#SBATCH --cpus-per-task=6
#SBATCH --job-name="threeloci_analysis"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out


###################################################################
#        Configuration 
###################################################################

# software requirements
module load StdEnv/2023 apptainer/1.3.5 seqkit/2.5.1

#  ENV variables for paths and input samplesheet
input_filepaths_list="/project/60006/mdprieto/giardia_mlst_2023/processed_data/accessions/fofn_HQclustered_classicMLST.txt"
mlst_fasta_dir="/project/60006/mdprieto/raw_data/giardia/MLST_loci/classic_three_loci"
results_dir="/scratch/mdprieto/results/giardia_classic_mlst/three_loci"

# ENV variables for containers and executable scripts
chewbbaca_img="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-chewbbaca-3.4.2-pyhdfd78af_0.sif"
iqtree_img="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-iqtree:3.0.1--h503566f_0.img"
seqkit_img="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-seqkit:2.10.1--he881be0_0.img"

###################################################################
#        Commands
###################################################################

# use the loci to create an adapted schema for chewBBACA
apptainer exec ${chewbbaca_img} chewBBACA.py PrepExternalSchema \
    --schema-directory "${mlst_fasta_dir}" \
    --translation-table 1 \
    --output-directory "${results_dir}/adapted_schema" \
    --cpu-cores 6

# perform allele calling using the adapted schema and extract the matching alleles as fasta
apptainer exec ${chewbbaca_img} chewBBACA.py AlleleCall \
    --input-files ${input_filepaths_list} \
    --schema-directory "${results_dir}/adapted_schema" \
    --translation-table 1 \
    --output-directory ${results_dir}/AlleleCall \
    --cpu-cores 6 &&
apptainer exec ${chewbbaca_img} chewBBACA.py GetAlleles \
    --input-file ${results_dir}/AlleleCall/results_alleles.tsv  \
    --schema-directory "${results_dir}/adapted_schema" \
    --translation-table 1 \
    --output-directory ${results_dir}/GetAlleles \
    --cpu-cores 6

# perform multiple sequence alignment for each locus
mkdir -p ${results_dir}/multiple_alignment &&
    for fasta_file in $(ls -1 ${results_dir}/GetAlleles/fastas/*.fasta)
    do  
        filename=$(basename -- "${fasta_file}")
        apptainer exec ${chewbbaca_img} mafft \
            --localpair \
            --maxiterate 1000 \
            "${fasta_file}" > "${results_dir}/multiple_alignment/${filename%.fasta}_aligned.fasta"
    done

# merge the multiple sequence alignments into a single concatenated alignment by sample
apptainer exec ${seqkit_img} seqkit concat \
    --full \
    --fill "-" \
    --id-regexp '(.*(?:SRR|ERR|SAMN|GCA)[0-9]+)' \
    --line-width 80 \
    --threads 6 \
    --out-file ${results_dir}/multiple_alignment/merged_MSA.fasta \
    ${results_dir}/multiple_alignment/NC*.fasta

# IQTREE2 with automatic model selection and 1000 ultrafast bootstraps + 1000 SH-aLRT replicates
mkdir -p "${results_dir}/iqtree" &&
    apptainer exec ${iqtree_img} iqtree2 \
        -s ${results_dir}/multiple_alignment/merged_MSA_three_loci.fasta \
        -m MFP \
        -B 1000 \
        -T AUTO \
        --alrt 1000  \
        --prefix "${results_dir}/iqtree/three_loci" \
        --verbose