#!/usr/bin/env bash

# ============================================================================= #
# PIPELINE CONFIGURATION                                                        #
# Centralized paths and environment variables for the poject.                   #
# Use with `source path/to/config_paths.sh`                                     #
# ============================================================================= #

# -------------------------------------------------------------------------
# EAGLE CONFIGURATION
# -------------------------------------------------------------------------

current_host=$(hostname)
drac_cluster="${CC_CLUSTER:-unknown}"

if [[ "$current_host" == *"eagle"* || "$current_host" =~ ^node[0-9]+ ]]; then

    # MANUAL: define Eagle-specific root paths
    project_root="/project/60006/mdprieto/giardia_mlst_2023"
    scratch_tmp_folder="/scratch/mdprieto/cache"

    # MANUAL: container path(s)
    sra_tools_container="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-sra-tools-3.2.1--h4304569_1.sif"

# -------------------------------------------------------------------------
# DRA CONFIGURATION
# -------------------------------------------------------------------------

elif [[ "$drac_cluster" == "fir" || "$current_host" == *"fir"* || "$current_host" =~ ^fc[0-9]+ ]]; then

    # MANUAL: define Fir-specific root paths
    project_root="/scratch/mdprieto/repositories/giardia_mlst_2023"
    scratch_tmp_folder="/scratch/mdprieto/cache"
    nf_chewbbacca_pipeline="/scratch/mdprieto/repositories/nf_chewbbaca_mlst"

    # MANUAL: container path(s) and environments
    sra_tools_container="/project/6007413/cidgoh_share/singularity_imgs/sra-tools-3.2.1--h4304569_1.img"
    quast_container="/project/6007413/cidgoh_share/singularity_imgs/quast-5.0.2--py37pl526hb5aa323_2.img"
    sourmash_container="/project/6007413/cidgoh_share/singularity_imgs/sourmash-4.8.9--hdfd78af_0.img"
    coverm_container="/project/6007413/cidgoh_share/singularity_imgs/coverm-0.7.0--hcb7b614_4.img"
    clustering_env="/home/mdprieto/virtual_envs/tsne_env"

    # MANUAL: databases
    eggnog_db_dir="/project/6007413/cidgoh_share/database/eagle/eggnog"

# FALLBACK/UNKNOWN
else
    echo "ERROR: Unknown computing environment. Host: ${current_host}, CC_CLUSTER: ${drac_cluster}"
    exit 1
fi

# -------------------------------------------------------------------------
# EXPORT ALL VARIABLES
# -------------------------------------------------------------------------

# Export manually defined path(s)
export project_root
export scratch_tmp_folder
export sra_tools_container
export quast_container
export sourmash_container
export coverm_container
export clustering_env
export eggnog_db_dir

# AUTOMATED: relative input to primary subfolders
export project_scripts_dir="${project_root}/scripts"
export input_dir="${project_root}/input_data"
export processed_data="${project_root}/processed_data"
export project_outdir="${project_root}/output"
export analysis_dir="${project_root}/analysis"

# AUTOMATED: relative path to other subfolders
export raw_reads_dir="${input_dir}/raw_fastq"
export reference_genome="${input_dir}/reference_genome"
export reference_genome_A="${reference_genome}/assemblage_A/GCF_000002435.2_WB_genomic.fna.gz"
export accessions_dir="${processed_data}/accessions"
export bactopia_results="${project_outdir}/bactopia_giardia_2025"
