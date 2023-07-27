#!/usr/bin/env nextflow

/*  
========================================================================================

    10 fold cross-validation of cgMLST calling in Giardia spp. data
========================================================================================
*/

//---------- Define params and input

params.folds=10
params.schema_path=""
params.outdir="$launchDir"

//---------- Processes


process CHEWBACCA_CREATE_SCHEMA {
    publishDir "${params.outdir}", mode: 'copy'
    cache 'lenient'
    cpus params.threads

    input:
        tuple val(sample_id), val(set), path(contig), val(train_test)
        path (prodigal_training)
    
    output:
        path("wgmlst_schema")

    script:
        """
        chewBBACA.py CreateSchema \
            -i $contigs.test -o $wgmlst_schema --ptf $traing_file --cpu $task.cpus
        """
}

process CHEWBACCA_ALLELE_CALL {
    publishDir "${params.outdir}", mode: 'copy'
    cache 'lenient'
    cpus params.threads

    input:
        tuple val(sample_id), val(set), path(contig), val(train_test)
        path (wgmlst_schema)
    
    output:
        path("results_AlleleCall")

    script:
        """
        chewBBACA.py AlleleCall \
            -i $contigs.test \
            -g $wgmlst_schema/schema_seed \
            -o results_AlleleCall \
            --cpu $task.cpus
        """
}

process CHEWBACCA_REMOVE_PARALOGS {
    publishDir "${params.outdir}", mode: 'copy'
    cache 'lenient'
    cpus params.threads

    input:
        path (results_AlleleCall)
    
    output:
        path("results_no_paralogs.tsv")

    script:
        """
        chewBBACA.py RemoveGenes \
            -i $results_AlleleCall/results_alleles.tsv \
            -g $results_AlleleCall/paralogous_counts.tsv \
            -o results_no_paralogs.tsv
        """
}

process CHEWBACCA_EXTRACT_CGMLST {
    publishDir "${params.outdir}", mode: 'copy'
    cache 'lenient'
    cpus params.threads

    input:
        tuple val(sample_id), val(set), path(contig), val(train_test)
        path (results_no_paralogs)
    
    output:
        path(cg)

    script:
        """
        chewBBACA.py RemoveGenes \
            -i $results_AlleleCall/results_alleles.tsv \
            -g $results_AlleleCall/paralogous_counts.tsv \
            -o results_no_paralogs.tsvq
        """
}
chewBBACA.py ExtractCgMLST \
    -i $results_no_paralogs \
    -o results32_wgMLST/cgMLST

//---------- Workflow