#!/usr/bin/env bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#                                                                               #
#       PIPELINE CONFIGURATION                                                  #
#           Centralized paths and environment variables for the poject.         #
#           Use with `source path/to/config_paths.sh`                           #
#                                                                               #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# =============================================================================
# Auto export and evaluate run environment
# =============================================================================

set -a

current_host=$(hostname)
drac_cluster="${CC_CLUSTER:-unknown}"

    # -------------------------------------------------------------------------
    # Eagle cluster (lab)
    # -------------------------------------------------------------------------

if [[ "$current_host" == *"eagle"* || "$current_host" =~ ^node[0-9]+ ]]; then

    # MANUAL: define Eagle-specific root paths
    project_root="/project/60006/mdprieto/giardia_mlst_2023"
    scratch_tmp_folder="/scratch/mdprieto/cache"

    # MANUAL: container path(s)
    sra_tools_container="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-sra-tools-3.2.1--h4304569_1.sif"

    # -------------------------------------------------------------------------
    # DRAC HPC
    # -------------------------------------------------------------------------

elif [[ "$drac_cluster" == "fir" || "$current_host" == *"fir"* || "$current_host" =~ ^fc[0-9]+ ]]; then

    # MANUAL: define Fir-specific root paths
    project_root="/scratch/mdprieto/repositories/giardia_mlst_2023"
    scratch_tmp_folder="/scratch/mdprieto/cache"
    nf_chewbbacca_pipeline="/scratch/mdprieto/repositories/nf_chewbbaca_mlst"
    nf_work_dir="/scratch/mdprieto/cache/nf_work_cache/"

    # MANUAL: databases
    eggnog_db_dir="/project/6007413/cidgoh_share/database/eagle/eggnog"
    kraken2_db="/project/6007413/cidgoh_share/database/eagle/Kraken2_bracken_db/kraken2_standard_20250402"

    # MANUAL: container path(s)
    coverm_container="/project/6007413/cidgoh_share/singularity_imgs/coverm-0.7.0--hcb7b614_4.img"
    sra_tools_container="/project/6007413/cidgoh_share/singularity_imgs/sra-tools-3.2.1--h4304569_1.img"
    kraken2_container="/project/6007413/cidgoh_share/singularity_imgs/community.wave.seqera.io-library-kraken2_coreutils_pigz-45764814c4bb5bf3.img"
    krakentools_container="/project/6007413/cidgoh_share/singularity_imgs/depot.galaxyproject.org-singularity-krakentools-1.2--pyh5e36f6f_0.img"
    quast_container="/project/6007413/cidgoh_share/singularity_imgs/quast-5.0.2--py37pl526hb5aa323_2.img"
    sourmash_container="/project/6007413/cidgoh_share/singularity_imgs/sourmash-4.8.9--hdfd78af_0.img"

    # MANUAL: python environment for clustering
    clustering_env="/home/mdprieto/virtual_envs/tsne_env2/"

    # -------------------------------------------------------------------------
    # Fail if unknown
    # -------------------------------------------------------------------------

else
    echo "ERROR: Unknown computing environment. Host: ${current_host}, CC_CLUSTER: ${drac_cluster}"
    exit 1
fi

# =============================================================================
# Automated sub-paths generation
# =============================================================================

# Major subfolders
project_scripts_dir="${project_root}/scripts"
input_dir="${project_root}/input_data"
processed_data="${project_root}/processed_data"
project_outdir="${project_root}/output"
analysis_dir="${project_root}/analysis"

# Other paths
raw_reads_dir="${input_dir}/raw_fastq"
reference_genome_dir="${input_dir}/reference_genome"
assemblage_A_fasta="${reference_genome_dir}/assemblage_A/GCF_000002435.2_WB_genomic.fna.gz"
accessions_dir="${processed_data}/accessions"
bactopia_results="${project_outdir}/bactopia_giardia_2025"

    # -------------------------------------------------------------------------
    # Kill auto-export
    # -------------------------------------------------------------------------

set +a
