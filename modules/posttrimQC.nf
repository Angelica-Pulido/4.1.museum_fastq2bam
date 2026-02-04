/*
   post-trimming_PE QC
*/

process posttrimQC_PE {
    label "fastqc"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC/2_posttrim", mode: params.publish

    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), path(paired), path(unpaired)

    output:
        tuple path("*.html"), path("*.zip")

    script:
    """
    module load ${params.fastqc}
    fastqc --noextract ${paired} ${unpaired}
    """

    stub:
    """
    touch ${basename}_stub.html
    touch ${basename}_stub.zip
    """
}

/*
   post-trimming_SE QC
*/

process posttrimQC_SE {
    label "fastqc"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC/2_posttrim", mode: params.publish

    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), path(assembledTrim)

    output:
        tuple path("*.html"), path("*.zip")

    script:
    """
    module load ${params.fastqc}
    fastqc --noextract ${assembledTrim}
    """

    stub:
    """
    touch ${basename}_stub.html
    touch ${basename}_stub.zip
    """
}

/*
   multiQC post-trimming
*/

process posttrim_multiQC {
    label "multiqc"
    tag "posttrim_multiqc"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC", mode: params.publish, pattern: "*html"
    publishDir "$params.result_dir/QC/2_posttrim", mode: params.publish, pattern: "*data"

    input:
        path(htmls)
        path(zips)

    output:
        file "**.html"
        file "*_data"

    script:
    """
    module load ${params.python}
    source ${params.python_venv}
    multiqc -n 2_posttrim_multiqc .
    """

    stub:
    """
    touch ${basename}.html
    touch ${basename}.zip
    """
}


workflow posttrim_qc_wf {
    take:
    channel_1
    channel_2

    main:
    posttrimQC_PE(channel_1)
    posttrimQC_SE(channel_2)

    // Collect all HTMLs and ZIPs from both processes
    all_htmls = posttrimQC_PE.out.mix(posttrimQC_SE.out).collect{ it[0] }
    all_zips = posttrimQC_PE.out.mix(posttrimQC_SE.out).collect{ it[1] }

    posttrim_multiQC(all_htmls, all_zips)
}