// Assign sex based on coverage

/*
    Samtools coverage 
*/

process coverage {
    label "coverage"
    tag "$source_organism_id"
    scratch params.scratch_dir
    publishDir "$params.result_dir/QC/6_coverage", mode: params.publish
    
    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), path(bam), path(bai)

    output:
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*.tsv")

    script:
    """
    dcsrsoft use ${params.softstack}
    module load ${params.samtools}

    samtools coverage ${bam} > ${specimen_id}_${lane}_coverage.tsv
    
    """
}

/*
   Sex assignment
*/

process sex_assign {
    label "sex_assign"
    tag "$source_organism_id"
    scratch params.scratch_dir
    // not necessary to publish per individual

    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), path(coverage)

    output:
        tuple val(source_organism_id), val(specimen_id), path("*.ssv")

    script:
    """
    dcsrsoft use ${params.softstack}
    module load ${params.python}
    source ${params.python_venv}
    
    assign_sex.py --source_organism_id ${source_organism_id} \
        --bam_id ${specimen_id} \
        --coverage ${coverage} \
        --ref-version ${params.ref_genome_version}
    """
}


workflow sex_assign_wf {
    take:
    markdupbams

    main:
    coverage(markdupbams)
    sex_ssv = sex_assign(coverage.out)
    // Collect all ssv into a single file
    sample_sex_assignment = sex_ssv
                            .map { it[2] }                   // grab the path only
                            .collectFile(name: 'sample_sex_assignment.txt', storeDir: "$params.result_dir")

    emit:
        sex_ssv
}
