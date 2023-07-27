#!/usr/bin/env nextflow

/*  
========================================================================================

    10 fold cross-validation of cgMLST calling in Giardia spp. data
========================================================================================
*/

// Input definition 

params.csv_files    = '/project/60006/mdprieto/giardia_mlst_2023/processed_data/cross_validation_input/pilot_input_long.csv'
params.ref_genome   = "/mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.fna"
params.organism     = "giardia_duodenalis"
params.outdir       = "$launchDir/results"

// Create channels

csv_channel = Channel.fromPath(params.csv_files)
ref_genome_ch = Channel.of([params.organism, params.ref_genome])

// Trial snippets

csv_channel
    .splitCsv(header:true)
    .map { return tuple(it.set, it.value, it.contig) }
    .set{input_channel}

def criteria = multiMapCriteria {
    set1_train: it[0] == "set1" && it[1] == "train"
    }

input_channel.multiMap(criteria).set{set1_train}
set1_train.view()