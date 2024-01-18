//------------- INPUT CHANNELS


params.csv_files = '../../processed_data/cross_validation_input/split_aa'
params.ref_genome = "/mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.fna"
params.organism = "giardia_duodenalis"

csv_channel = Channel.fromPath(params.csv_files)
ref_genome_ch = Channel.of([params.organism, params.ref_genome])


//------------- PROCESSES

process PRODIGAL_TRAINING {
    label "process_single"
    container "https://depot.galaxyproject.org/singularity/prodigal%3A2.6.3--hec16e2b_5"


    input:
    tuple val(organism), path (ref_genome)

    output:
    tuple val(organism), path ("*.trn"),    emit: training_file

    script:
    """
    prodigal \
        -i $ref_genome \
        -t giardia_wb.trn \
        -p single
    """
}

process TEST_CREATE_SCHEMA {
    //publishDir "${params.outdir}", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"
    debug true

    input:
    path(contigs_test)

    script:
    """
    realpath $contigs_test > contigs_input.txt
    cat contigs_input.txt
    """
}

process CHEWBACCA_CREATE_SCHEMA {
    //publishDir "${params.outdir}", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"

    input:
    path(contigs_test)
    tuple val(organism), path (training_file)
    
    output:
    path("wgmlst_schema")

    script:
    """
    realpath $contigs_test > contigs_input.txt
    cat contigs_input.txt


    chewBBACA.py CreateSchema \
        -i $contigs_test -o wgmlst_schema --ptf $training_file --cpu $task.cpus
    """
}

workflow{
    // create contigs_channel
    csv_channel
        .splitCsv(header:true)
        .branch { 
            train: it.value == "train"
                return it.contig
            test: it.value == "test"
                return it.contig  }
        .set{contigs_channel}

    training_ch = PRODIGAL_TRAINING(ref_genome_ch)
    
    TEST_CREATE_SCHEMA(contigs_channel.test.flatten().collect())
    //CHEWBACCA_CREATE_SCHEMA(contigs_channel.test.collect(), training_ch)

}

/*
//-------- SNIPPETS TO TRY

process CREATE_SET_CHANNELS{

    input:
    path (csv_file)

    output:
    stdout

    script:
    """
    #!/usr/bin/env groovy

    $csv_file
        .splitCsv(header:true)
        .branch { row ->
            train: row.value == "train"
                return tuple(row.sample, row.set, row.contig)
            test: row.value == "test"
                return tuple(row.sample, row.set, row.contig)
                }
        .set{set_ch}
    
    set_ch.train.view()
    """
}


*/