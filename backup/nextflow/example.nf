    // take all data and branch in subsets train/test
Channel.fromPath('example.csv')
    .splitCsv(header:true)
    .branch {   train: it.value == "train"
                    return tuple(it.set, it.contig)
                test: it.value == "test"
                    return tuple(it.set, it.contig)  }
    .set{contigs_channel}

process FOO {
    debug true

    input:
    tuple val(set), path(contigs_test)

    script:
    """
    echo $contigs_test | tr -s ' ' '\n' > contigs_input.txt
    cat contigs_input.txt
    
    echo command -i contigs_input.txt -o output_dir -t training_file
    """
}

workflow{
        // group it by set_id
    train_contigs_ch = contigs_channel.train.groupTuple()

    FOO(train_contigs_ch)
}

