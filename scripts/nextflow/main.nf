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