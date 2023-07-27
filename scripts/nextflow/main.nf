/*

// -----------  Trials of input definition


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

// -----------  Alternative contigs_channel

    csv_channel
        .splitCsv(header:true)
        .filter {it.value == "test"}        // filter out training subset
        .map {return it.contig }            // take only the PATH to contigs
        .collect().flatten()                // make it so every file is in a different line
        .take(-1)
        .set {test_contigs_ch}
    
    test_contigs_ch.view()

*/

// Input definition 

params.csv_files    = '../../processed_data/cross_validation_input/split_aa'
params.ref_genome   = "/mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.fna"
params.organism     = "giardia_duodenalis"

// Channel creation 

csv_channel = Channel.fromPath(params.csv_files)
ref_genome_ch = Channel.of([params.organism, params.ref_genome])

// Processes

process PRODIGAL_TRAINING {
    label "process_single"
    container "https://depot.galaxyproject.org/singularity/prodigal%3A2.6.3--hec16e2b_5"


    input:
    tuple val(organism), path (ref_genome)

    output:
    tuple val(organism), path ("*.trn")

    script:
    """
    prodigal \
        -i $ref_genome \
        -t giardia_wb.trn \
        -p single
    """
}

process TEST_PROCESS {
    debug true

    input:
    tuple val(set_id), path (contigs_test)

    output:
    stdout

    script:
    """
    echo $contigs_test | tr -s ' ' '\n' > contigs_test.txt
    cat contigs_test.txt
    """

}

process CHEWBACCA_CREATE_SCHEMA {
    //publishDir "${params.outdir}", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"
    cpus 6

    input:
    tuple val(set_id), path (contigs_test)
    tuple val(organism), path (training_file)
    
    output:
    path("wgmlst_schema")

    script:
    """

    echo $contigs_test | tr -s ' ' '\n' > contigs_input.txt
    cat contigs_input.txt

    chewBBACA.py CreateSchema \
        -i contigs_input.txt -o wgmlst_schema --ptf $training_file --cpu $task.cpus
    """
}

workflow{
    csv_channel
        .splitCsv(header:true)
        .branch { 
            train: it.value == "train"
                return tuple(it.set, it.contig)
            test: it.value == "test"
                return tuple(it.set, it.contig)  }
        .set{contigs_channel}

    contigs_channel.train
        .take(20)
        .groupTuple()
        .set{ train_channel }

    training_ch = PRODIGAL_TRAINING(ref_genome_ch)
    // TEST_PROCESS(train_channel)
    CHEWBACCA_CREATE_SCHEMA(
        train_channel,
        training_ch
    )
    }

    // contigs_channel.train.view()




