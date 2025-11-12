#!/bin/bash
#SBATCH --mem-per-cpu=3G
#SBATCH --time=3-05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --job-name="bactopia_hybrid"
#SBATCH --chdir=/scratch/mdprieto/
#SBATCH --output=jobs_output/%j_%x.out

###############################################################################################

# load modules/dependencies
module load nextflow apptainer

# define environment variables for HPC
project_root="/project/60006/mdprieto/giardia_mlst_2023"
hydrid_fastq_dir="${project_root}/input_data/giardia/repositories/fastq"
custom_config="${project_root}scripts/data_retrieval_assembly/eagle_bactopia.config"

###############################################################################################

    # no samplesheet option available for hybrid runs in Bactopia 3.0.0

nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    --sample SAMN12610744 \
    --r1 ${hydrid_fastq_dir}/SRX6745910_SRR10007607_1.fastq.gz \
    --r2 ${hydrid_fastq_dir}/SRX6745910_SRR10007607_2.fastq.gz \
    --ont ${hydrid_fastq_dir}/SRX6745908-SRR10007609.fastq.gz \
    -resume \
    --nfconfig ${custom_config} \
    --outdir /scratch/mdprieto/results/bactopia_giardia/hybrid/ \
    --short_polish \
    --skip_amr \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 180040666

nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    --sample SAMN12611509 \
    --r1 ${hydrid_fastq_dir}/SRX6745954_SRR10007650_1.fastq.gz \
    --r2 ${hydrid_fastq_dir}/SRX6745954_SRR10007650_2.fastq.gz \
    --ont ${hydrid_fastq_dir}/SRX6745952-SRR10007652.fastq.gz \
    -resume \
    --nfconfig ${custom_config} \
    --outdir /scratch/mdprieto/results/bactopia_giardia/hybrid/  \
    --short_polish \
    --skip_amr \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 180040666

nextflow run bactopia/bactopia -r v3.0.0 \
    -profile singularity,slurm \
    --sample SAMN12611599 \
    --r1 ${hydrid_fastq_dir}/SRX6746024_SRR10007722_1.fastq.gz \
    --r2 ${hydrid_fastq_dir}/SRX6746024_SRR10007722_2.fastq.gz \
    --ont ${hydrid_fastq_dir}/SRX6746022-SRR10007724.fastq.gz \
    -resume \
    --nfconfig ${custom_config} \
    --outdir /scratch/mdprieto/results/bactopia_giardia/hybrid/  \
    --short_polish \
    --skip_amr \
    --skip-prokka \
    --skip_mlst \
    --cleanup_workdir \
    --max_genome_size 180040666
