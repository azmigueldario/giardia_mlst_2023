# requires apptainer

# define PATH with fastq files
FASTQ_DIR='/project/60006/mdprieto/raw_data/giardia/repositories/fastq'

# remove low hyphen from ONT reads (necessary before bactopia prepare)
for file in $FASTQ_DIR/*fastq.gz
    do
    echo mv $file $(echo $file | sed -E 's|([0-9]*)_([SRR,ERR])|\1-\2|')
    done

# requires Bactopia apptainer image from https://depot.galaxyproject.org/singularity/

# run bactopia prepare script to produce samplesheet
singularity exec -B /etc /mnt/cidgoh-object-storage/images/bactopia_v3.0.0.sif \
    bactopia prepare \
    --path $FASTQ_DIR \
    --ont > bactopia_samplesheeet_eagle.csv