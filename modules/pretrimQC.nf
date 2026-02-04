// Pre-trimming QC subworkflow

/*
   FastQC pre-trimming
*/

process pretrimQC {
    label "fastqc"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC/0_pretrim", mode: params.publish
    
    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), file(r1), file(r2)

    output:
        tuple path("*.html"), path("*.zip")

    script:
    """
    module load ${params.fastqc}
    fastqc --noextract ${r1} ${r2} 
    """
}

/*
   multiQC pre-trimming
*/

process pretrim_multiQC {
    label "multiqc"
    tag "pretrim_multiqc"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC", mode: params.publish, pattern: "*html"
    publishDir "$params.result_dir/QC/0_pretrim", mode: params.publish, pattern: "*data"

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
    multiqc -n 0_pretrim_multiqc .
    """

    stub:
    """
    touch ${basename}.html
    touch ${basename}.zip
    """
}

workflow pretrim_qc_wf {
    take:
    renamed_reads

    main:
    pretrimQC(renamed_reads)
    
    // Collect and pass to multiqc
    all_htmls = pretrimQC.out.collect{ it[0] }  // Collect all HTML files
    all_zips = pretrimQC.out.collect{ it[1] }   // Collect all ZIP files
    
    pretrim_multiQC(all_htmls, all_zips)
}