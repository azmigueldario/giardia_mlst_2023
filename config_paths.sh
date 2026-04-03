#!/usr/bin/env bash

# ==============================================================================
# PIPELINE CONFIGURATION
# Description: Centralized paths and environment variables for the workflow.
# Usage: source path/to/config.sh
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Base paths [Manually defined]
# ------------------------------------------------------------------------------
project_root="/project/60006/mdprieto/giardia_mlst_2023"
scratch_tmp_folder="/scratch/mdprieto/cache"
project_scripts_dir="${project_root}/scripts"
processed_data="${project_root}/processed_data"

# ------------------------------------------------------------------------------
# 2. Container path(s) [Manually defined]
# ------------------------------------------------------------------------------
sra_tools_container="/scratch/group_share/singularity_imgs/depot.galaxyproject.org-singularity-sra-tools-3.2.1--h4304569_1.sif"

# ------------------------------------------------------------------------------
# 3. Input data paths
# ------------------------------------------------------------------------------
raw_reads_dir="${project_root}/input_data/raw_fastq"
reference_genome="${project_root}/input_data/reference_assemblage_A"
accessions_dir="${processed_data}/accessions"

# ------------------------------------------------------------------------------
# 4. Output destinations
# ------------------------------------------------------------------------------
project_outdir="${project_root}/output"
