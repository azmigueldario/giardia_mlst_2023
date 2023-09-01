 # to create the bactopia input samplesheet
OUTPUT_DIR="PATH/TO/OUTPUT/DIR"

 echo -e "sample\truntype\tr1\tr2\textra" > $OUTPUT_DIR/bactopia_samplesheet.csv

for read1 in $(ls /project/6056895/mdprieto/raw_data/giardia/{BCCDC,repositories}/fastq/*_1.fastq.gz);
    do
        # get the basename of file and remove suffix
    sample_id=$(echo $read1 | xargs -n 1 basename -s '_1.fastq.gz')
        # replace string '_1' for '_2'
    read2="${read1/_1/_2}"
        # write in a new line for each sample
    echo -e "$sample_id\tpaired-end\t$read1\t$read2"  >> $OUTPUT_DIR/bactopia_samplesheet.csv
    done 
