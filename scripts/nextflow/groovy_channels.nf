// -----------  Groovy function to create channels
def create_channels(LinkedHashMap row) {
    
    def meta = [:]              // create meta map
    meta.id         = row.sample

    def fastq_meta = []
    if (mode.matches('nanopore')) {
        fastq_meta = [ meta, [ file(row.fastq_dir) ] ]
    }
    else {
        if (!file(row.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        if (meta.single_end) {
            fastq_meta = [ meta, [ file(row.fastq_1) ] ]
        } else {
            if (!file(row.fastq_2).exists()) {
                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
            }
            fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
        }
    }
    return fastq_meta
}Z