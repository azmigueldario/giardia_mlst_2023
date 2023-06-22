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


process SUBSET_CROSS_VALIDATION {
    container 
    label 'low'
    publishDir "${params.outdir}/", mode: 'symlink'
}

process CHEW:ALLELE_CALL{
    tag "Processing $sample_id"
    container 
    label 'low'
    publishDir "${params.outdir}/{}", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id)
    


}

//---------- Workflow