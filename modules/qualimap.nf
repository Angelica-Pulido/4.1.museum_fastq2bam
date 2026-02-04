// Qualimap qc workflow

/*
   Qualimap
*/

process qualimap {
    label "qualimap"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC/4_bamqc", mode: params.publish
    
    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), path(bam), path(bai)
        path bed

    output:
        file "*_bamqc"

    script:
    """
    module load ${params.gcc} ${params.openjdk} ${params.qualimap}
    new_name=${specimen_id}_${params.run_id}_${lane}
    export JAVA_OPTS="-Xmx${task.memory.giga}G"

    qualimap bamqc -c -bam ${bam} -gff ${bed} -nw ${params.qualimap_windows} -outdir \${new_name}_bamqc -nt ${task.cpus} -outformat PDF:HTML
    """
}

/*
   multiQC pre-trimming
*/

process qualimap_multiQC {
    label "multiqc"
    tag "qualimap_multiqc"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/QC", mode: params.publish, pattern: "*html"
    publishDir "$params.result_dir/QC/4_bamqc", mode: params.publish, pattern: "*data"

    input:
        path(bamqcdir)

    output:
        file "**.html"
        file "*_data"

    script:
    """
    module load ${params.python}
    source ${params.python_venv}
    multiqc -n 4_bamqc_multiqc .
    """
}

workflow qualimap_wf {
    take:
    bam
    bed

    main:
    qualimap(bam, bed) 
    qualimap_multiQC(qualimap.out.collect())
}