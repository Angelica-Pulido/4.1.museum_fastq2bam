/*
   multiQC of trimmomatic logs
*/

process trim_stats {
    label "multiqc"
    tag "trim_multiqc"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC", mode: params.publish, pattern: "*html"
    publishDir "$params.result_dir/QC/1_trimstats", mode: params.publish, pattern: "*data"

    input:
        path(logs)

    output:
        file "**.html"
        file "*_data"

    script:
    """
    module load ${params.python}
    source ${params.python_venv}
    multiqc -n trim_multiqc .
    """

    stub:
    """
    touch ${basename}.html
    touch ${basename}.zip
    """
}