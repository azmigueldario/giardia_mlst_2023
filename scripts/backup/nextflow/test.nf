#!/usr/bin/env nextflow

/*  
========================================================================================

    10 fold cross-validation of cgMLST calling in Giardia spp. data
========================================================================================
*/


// Input definition 

    // csv_input may have one subset or multiple, here it has three
params.csv_files    = '/project/60006/mdprieto/giardia_mlst_2023/processed_data/cross_validation_input/pilot_input_long.csv'
params.ref_genome   = "/mnt/cidgoh-object-storage/database/reference_genomes/giardia/assemblage_A/GCF_000002435.2_UU_WB_2.1_genomic.fna"
params.organism     = "giardia_duodenalis"
params.outdir       = "$launchDir/results"

// Create channels

csv_channel = Channel.fromPath(params.csv_files)
ref_genome_ch = Channel.of([params.organism, params.ref_genome])

// --------------- Trial snippets

// Branch criteria: define criteria for separation and apply to input channel after splitCsv

def criteria_branch = branchCriteria {it ->        
    train: it.value=="train"
        return tuple(it.set, it.contig)
    test: it.value == "test"
        return tuple(it.set, it.contig)  }

    csv_channel
        .splitCsv(header:true)
        .branch(criteria_branch)
        .set{contigs_channel}

// Use groupTuple() to make a tuple containing the identifier value and a list of all contigs of a subset
    training_channel = contigs_channel.train.groupTuple()
    testing_channel  = contigs_channel.test.groupTuple()
    
process FOO{
    debug true
    
    input:
    tuple val(set_id), path (training_contigs)
    
    script:
    """
        # creates a file with one contig absolute path per line
    echo $training_contigs | tr -s ' ' '\n' > train_contigs_list.txt
    
    cat train_contigs_list.txt
    """
}

workflow{
    FOO(training_channel)
    }
/* GroupTuple directly
    csv_channel
        .splitCsv(header:true)
        .map(row -> tuple(row.set, row.value, row.contig) )
        .groupTuple(by: [0,1])
        .set{contigs_channel}
    
    contigs_channel.view()


/*  
========================================================================================

                    Working as expected
========================================================================================

// Branch criteria: define criteria for separation and apply to input channel after splitCsv

def criteria_branch = branchCriteria {it ->        
    train: it.value=="train"
        return tuple(it.set, it.contig)
    test: it.value == "test"
        return tuple(it.set, it.contig)  }

    csv_channel
        .splitCsv(header:true)
        .branch(criteria_branch)
        .set{contigs_channel}

// Use groupTuple() to make a tuple containing the identifier value and a list of all contigs of a subset
    training_channel = contigs_channel.train.groupTuple()
    testing_channel  = contigs_channel.test.groupTuple()

// Use join() to mix together the contigs belonging to a same set.
// Joining fails if duplicates are found or if there are no matches in the set.
    training_channel
        .join(testing_channel, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .view()

========================================================================================

                    Not working yet
========================================================================================


    // Multimap
csv_channel
    .splitCsv(header:true)
    .map { return tuple(it.set, it.value, it.contig) }
    .set{input_channel}

def criteria = multiMapCriteria { it ->
    set1_train: it[0] == "set1" && it[1] == "train"
    other: true
    }
input_channel.multiMap(criteria).set{output_ch}
output_ch.set1_train.view{it}
*/