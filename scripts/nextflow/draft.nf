#!/usr/bin/env nextflow

/*  
========================================================================================

    10 fold cross-validation of cgMLST calling in Giardia spp. data
========================================================================================
*/

// Input definition 

params.csv_files    = '../../processed_data/cross_validation_input/split_aa'
params.ref_genome   = "/mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.fna"
params.organism     = "giardia_duodenalis"
params.outdir       = "/scratch/mdprieto/results/giardia_chewBBACA"//"$launchDir/results"

// Create channels

csv_channel = Channel.fromPath(params.csv_files)
ref_genome_ch = Channel.of([params.organism, params.ref_genome])

// Processes

process PRODIGAL_TRAINING {
    label "process_single"
    container "https://depot.galaxyproject.org/singularity/prodigal%3A2.6.3--hec16e2b_5"
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    tuple val(organism), path (ref_genome)

    output:
    tuple val(organism), path ("*.trn")

    script:
    """
    prodigal \
        -i $ref_genome \
        -t ${organism}.trn \
        -p single
    """
}

process CHEWBACCA_CREATE_SCHEMA {
    //label "process_medium"
    cpus 8
    publishDir "${params.outdir}/${set_id}", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"

    input:
    tuple val(set_id), path (training_contigs)
    tuple val(organism), path (prodigal_training)
    
    output:
    tuple val(set_id), path("wgmlst_schema"), path('train_contigs_list.txt')
    
    script:
    """
        # creates a file with one contig absolute path per line
    echo $training_contigs | tr -s ' ' '\n' > train_contigs_list.txt

    chewBBACA.py CreateSchema \
        -i train_contigs_list.txt \
        -o wgmlst_schema \
        --ptf $prodigal_training \
        --cpu $task.cpus
    """
}

process CHEWBACCA_ALLELE_CALL {
    label 'process_medium'
    publishDir "${params.outdir}/${set_id}/train/", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"

    input:
    tuple val(set_id), path(wgmlst_schema), path(train_contigs_list)
    tuple val(set_id), path (training_contigs)
    
    output:
    tuple val(set_id), path("results_AlleleCall")

    script:
    """
    chewBBACA.py AlleleCall \
        -i $train_contigs_list \
        -g $wgmlst_schema/schema_seed \
        -o results_AlleleCall \
        --cpu $task.cpus
    """
}

process REMOVE_PARALOGS {
    label 'process_medium'
    publishDir "${params.outdir}/${set_id}", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"

    input:
    tuple val(set_id), path(results_AlleleCall)
    
    output:
    tuple val(set_id), path("**results_no_paralogs.tsv")

    script:
    """
    chewBBACA.py RemoveGenes \
        -i $results_AlleleCall/results_alleles.tsv \
        -g $results_AlleleCalln /paralogous_counts.tsv \
        -o ${set_id}_results_no_paralogs_train.tsv
    """
}

process EXTRACT_CGMLST {
    label 'process_medium'
    publishDir "${params.outdir}/${set_id}", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"

    input:
    tuple val(set_id), path(no_paralogs_results)
    
    output:
    tuple val(set_id), path('*cgMLST')

    script:
    """
    chewBBACA.py ExtractCgMLST \
        -i $no_paralogs_results \
        -o train_cgMLST
    """
}

process ALLELE_CALL_TEST {
    label 'process_medium'
    publishDir "${params.outdir}/${set_id}", mode: 'copy'
    cache 'lenient'
    container "https://depot.galaxyproject.org/singularity/chewbbaca%3A3.2.0--pyhdfd78af_0"

    input:
    tuple val(set_id), path(wgmlst_schema), path(test_contigs)
    tuple val(set_id), path (contigs_train)
    
    output:
    tuple val(set_id), path("results_AlleleCall")

    script:
    """
        # creates a file with one contig absolute path per line
    echo $contigs_test | tr -s ' ' '\n' > test_contigs_input.txt

    chewBBACA.py AlleleCall \
        -i $test_contigs \
        -g $wgmlst_schema/schema_seed \
        -o results_AlleleCall \
        --cpu $task.cpus
    """
}

workflow{
    csv_channel
        .splitCsv(header:true)
        .branch { train: it.value == "train"
                    return tuple(it.set, it.contig)
                  test: it.value == "test"
                    return tuple(it.set, it.contig)  }
        .set{contigs_channel}

    contigs_channel.train
        .take(20)
        .groupTuple().set{ train_channel }
    
    contigs_channel.test
        .take(5)
        .groupTuple().set{ test_channel }
    
    prodigal_ch     = PRODIGAL_TRAINING(ref_genome_ch)
    wgmlst_ch       = CHEWBACCA_CREATE_SCHEMA(train_channel, prodigal_ch)
    wgmlst_schema   = CHEWBACCA_ALLELE_CALL(wgmlst_ch, train_channel)
    REMOVE_PARALOGS(wgmlst_schema)
    EXTRACT_CGMLST(REMOVE_PARALOGS.out)

    }


/*

//---------- Define params and input

params.folds=10
params.schema_path=""
params.outdir="$launchDir"

// ------------- Mix training and testing datasets once again
train_channel
    .join(test_channel, by: 0, failOnDuplicate: true, failOnMismatch: true)
    .view()

*/
