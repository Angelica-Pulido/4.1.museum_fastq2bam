/*
   multiQC marked dups metrics
*/

process dups_multiQC {
    label "multiqc"
    tag "dups_multiqc"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC", mode: params.publish, pattern: "*html"
    publishDir "$params.result_dir/QC/5_markedups", mode: params.publish, pattern: "*data"

    input:
        path(bamqcdir)

    output:
        file "**.html"
        file "*_data"

    script:
    """
    module load ${params.python}
    source ${params.python_venv}
    multiqc -n 5_dups_multiqc .
    """

    stub:
    """
    touch ${basename}.html
    touch ${basename}.zip
    """
}