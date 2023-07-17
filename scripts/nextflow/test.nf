/*

// Trying to do it natively through branch based on set name

Channel
    .fromPath(params.input_file)
    .splitCsv(header: true)
    .branch{ row ->
        set1: row.set1_train == "train"
            return tuple(row.sample, "set1", tuple(row.read1, row.read2))
        set2: row.set1_test == "test"
            return tuple(row.sample, "set1", tuple(row.read1, row.read2))  
    }
    .set{test_ch}

test_ch.set1.view{it}
test_ch.set2.view{it}


// input every file individually and run the pipeline 10 times

Channel
    .fromPath('path/to/set1_train.csv')
    .splitCsv(header: true)
    .map {row -> tuple(row.sample_id, row.set, row.contig) }
    .set {train_files_ch}

Channel
    .fromPath('path/to/set1_test.csv')
    .splitCsv(header: true)
    .map {row -> tuple(row.sample_id, row.set, row.contig) }
    .set {test_files_ch}

// Groovy function to create channels
def create_channels(LinkedHashMap row) {
    
    def meta = [:]              // create meta map
    meta.id         = row.sample

    // add path(s) of the fastq file(s) to the meta map
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
}

*/

params.input_file = "../../processed_data/cross_validation_input/pilot_cv.csv"

Channel
    .fromPath(params.input_file)
    .splitCsv(header: true)
    .branch { line -> 
            set1: 
            tuple(row.sample, "set1", row.set1, row.contig) 
            set2: tuple}
    .view()