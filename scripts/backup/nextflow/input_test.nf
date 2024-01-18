
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
    tuple val(organism), path ("*.trn")

    script:
    """
    prodigal \
        -i $ref_genome \
        -t giardia_wb.trn \
        -p single
    """

}

/*
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

workflow{
    
    training_ch = PRODIGAL_TRAINING(ref_genome_ch)

    csv_channel
        .splitCsv(header:true)
        .branch {train: it.value == "train"
                     return tuple(it.sample, it.set, it.contig)
                 test: it.value == "test"
                     return tuple(it.sample, it.set, it.contig)}
        .set{set_ch}

    set_ch.train
        .map{it[2]}
        .take (40)
        .collect().flatten()

    //test groupTuple
    csv_channel
        .splitCsv(header:true)
        .map {return tuple(it.set, it.contig)}
        .groupTuple()
        .view()
    //set_ch.train.groupTuple{}

    //ref_genome_ch.view(); csv_channel.view(); training_ch.view()
}