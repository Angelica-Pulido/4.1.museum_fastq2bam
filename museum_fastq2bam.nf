/*
   Get read groups
*/

process getRG {
    label "getRG"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    // this does not need to be published; uncomment otherwise
    // publishDir "$params.result_dir/RGgroups", mode: params.publish

    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), file(r1), file(r2)

    output:
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*.txt")

    script:
    """
    dcsrsoft use ${params.softstack}
    module load ${params.python}
    source ${params.python_venv}
    ${params.bin_dir}/rg_parser.py --source_organism_id ${source_organism_id} --specimen_id ${specimen_id} --run_id ${params.run_id} --lane ${lane} --r1 ${r1}
    """
}

/*
    Merge overlaping reads
    FASTQs will be renamed following the convention: {specimen_id}_{run_id}_{lane}
    */
    //NOTE: this process run and generared the files but didn't publish them!!! it's because it had: publishDir "$params.result_dir/1_trimmed", mode: params.publish, pattern: "*fastq.gz"

process merge_overlaping_reads {
    label "merge"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/1_trimmed", mode: params.publish, pattern: "*assembled.fastq"

    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), file(r1), file(r2)

    output:
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*assembled.fastq")

    script:
    """
    dcsrsoft use ${params.softstack}
    new_name=${specimen_id}_${params.run_id}_${lane}

    ${params.bin_dir}/${params.pear} \\
    -j ${task.cpus} \\
    -f ${r1} \\
    -r ${r2} \\
    -o \${new_name}
       
    """
}

/*
   Trim fastq
   FASTQs will be renamed following the convention: {specimen_id}_{run_id}_{lane}
*/
// NOTE: what is this stub?? at the end for?

process trim_PE {
    label "trim_PE"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/1_trimmed", mode: params.publish, pattern: "*fastq.gz"
    publishDir "$params.result_dir/QC/1_trimstats", mode: params.publish, pattern: "*.log"

    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), file(r1), file(r2)
        path adapter_file

    output:
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*_paired.fastq.gz"), path("*_unpaired.fastq.gz"), emit: trimmed_PE
        path("*trimPE.log"), emit: trimstats

    script:
    """
    dcsrsoft use ${params.softstack}
    module load  ${params.trimmomatic}
    new_name=${specimen_id}_${params.run_id}_${lane}

    trimmomatic PE \\
        -threads ${task.cpus} \\
        -${params.phred} \\
        ${r1} \\
        ${r2} \\
        \${new_name}_R1_paired.fastq.gz \\
        \${new_name}_R1_unpaired.fastq.gz \\
        \${new_name}_R2_paired.fastq.gz \\
        \${new_name}_R2_unpaired.fastq.gz \\
        ILLUMINACLIP:${adapter_file}:${params.illumina_clipping} \\
    	MINLEN:${params.min_read_length} 2> \${new_name}_trimPE.log 
    """
}
    
process trim_SE {
    label "trim_SE"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/1_trimmed", mode: params.publish, pattern: "*fastq.gz"
    publishDir "$params.result_dir/QC/1_trimstats", mode: params.publish, pattern: "*.log"

    input:
        tuple val(source_organism_id), val(specimen_id), val(lane), file(assembled)
        path adapter_file

    output:
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*_assembledTrim.fastq.gz"), emit: trimmed_SE
        path("*trimSE.log"), emit: trimstats

    script:
    """
    dcsrsoft use ${params.softstack}
    module load  ${params.trimmomatic}
    new_name=${specimen_id}_${params.run_id}_${lane}
    
    trimmomatic SE \\
        -threads ${task.cpus} \\
        -${params.phred} \\
        ${assembled} \\
        \${new_name}_assembledTrim.fastq.gz \\
        ILLUMINACLIP:${adapter_file}:${params.illumina_clipping} \\
        MINLEN:${params.min_read_length} 2> \${new_name}_trimSE.log
        
    """
}
/*
    stub:
    """
    new_name=\$(basename ${reads[0]} | sed 's/_[12]\\.fastq\\.gz//')
    touch \${new_name}_paired.fastq.gz
    touch \${new_name}_unpaired.fastq.gz
    touch \${new_name}_assembledTrim.fq.gz
    touch \${new_name}_trim.log
    """
/*
/*
    Map to reference genome
*/

process bwa_map {
    label "bwa_map"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    // this does not need to be published; uncomment otherwise
    // publishDir "$params.result_dir/2_mapped", mode: params.publish, pattern: "*.sam"
    publishDir "$params.result_dir/QC/3_bwalogs", mode: params.publish, pattern: "*_bwa.log"

    input: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path(paired_reads), path(unpaired_reads), path(assembled_reads)
        path ref_index

    output: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*_PE_mapped.sam"), path("*_UPF_mapped.sam"), path("*_UPR_mapped.sam"), path("*_SE_mapped.sam"), emit: aln
        path("*_bwa.log"), emit: log   
        // test and report the type of output
    script: 
        def idxbase = ref_index[0].baseName
        // NOTE!!!: why not putting the ref fasta file instead of the index base?
    """
    dcsrsoft use ${params.softstack}
    module load ${params.bwa}
    new_name=${specimen_id}_${params.run_id}_${lane}

    bwa mem -M -t ${task.cpus} ${idxbase} ${paired_reads[0]} ${paired_reads[1]} > \${new_name}_PE_mapped.sam 2>  \${new_name}_PE_mapped_bwa.log
    bwa mem -M -t ${task.cpus} ${idxbase} ${unpaired_reads[0]} > \${new_name}_UPF_mapped.sam 2>  \${new_name}_UPF_mapped_bwa.log
    bwa mem -M -t ${task.cpus} ${idxbase} ${unpaired_reads[1]} > \${new_name}_UPR_mapped.sam 2>  \${new_name}_UPR_mapped_bwa.log
    bwa mem -M -t ${task.cpus} ${idxbase} ${assembled_reads} > \${new_name}_SE_mapped.sam 2>  \${new_name}_SE_mapped_bwa.log
    """
}

/*
    Add read groups
*/

process read_groups {
    label "readgroups"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    // this does not need to be published; uncomment otherwise
    // publishDir "$params.result_dir/2_mapped", mode: params.publish

    input: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path(PE), path(UPF), path(UPR), path(SE), path(rg_info)
    output: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*_PE_RG.sam"), path("*_UPF_RG.sam"), path("*_UPR_RG.sam"), path("*_SE_RG.sam")

    script: 
    
    """
    dcsrsoft use ${params.softstack}
    module load ${params.samtools}

    new_name=${specimen_id}_${params.run_id}_${lane}

    lib=\$(cat ${rg_info} | awk -F"|" '{print \$1}')
    ID=\$(cat ${rg_info} | awk -F"|" '{print \$2}')
    BC=\$(cat ${rg_info} | awk -F"|" '{print \$3}')
    SM=\$(cat ${rg_info} | awk -F"|" '{print \$4}')

    samtools addreplacerg -r "ID:\$ID" \
        -r "BC:\$BC" \
        -r "CN:${params.RG_CN}" \
        -r "LB:\$lib" \
        -r "PL:${params.RG_PL}" \
        -r "PM:${params.RG_PM}" \
        -r "SM:\$SM" \
        -@ ${task.cpus} \
        -o \${new_name}_PE_RG.sam \
        ${PE}

    samtools addreplacerg -r "ID:\$ID" \
        -r "BC:\$BC" \
        -r "CN:${params.RG_CN}" \
        -r "LB:\$lib" \
        -r "PL:${params.RG_PL}" \
        -r "PM:${params.RG_PM}" \
        -r "SM:\$SM" \
        -@ ${task.cpus} \
        -o \${new_name}_UPF_RG.sam \
        ${UPF}

    samtools addreplacerg -r "ID:\$ID" \
        -r "BC:\$BC" \
        -r "CN:${params.RG_CN}" \
        -r "LB:\$lib" \
        -r "PL:${params.RG_PL}" \
        -r "PM:${params.RG_PM}" \
        -r "SM:\$SM" \
        -@ ${task.cpus} \
        -o \${new_name}_UPR_RG.sam \
        ${UPR}

    samtools addreplacerg -r "ID:\$ID" \
        -r "BC:\$BC" \
        -r "CN:${params.RG_CN}" \
        -r "LB:\$lib" \
        -r "PL:${params.RG_PL}" \
        -r "PM:${params.RG_PM}" \
        -r "SM:\$SM" \
        -@ ${task.cpus} \
        -o \${new_name}_SE_RG.sam \
        ${SE}

    """
}


/*
    Sort reads & output BAM
*/

process sortbam {
    label "sortbam"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    // this does not need to be published; uncomment otherwise
    // publishDir "$params.result_dir/2_mapped", mode: params.publish

    input: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path(PE), path(UPF), path(UPR), path(SE)
    output: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*_PE_sorted.bam"), path("*_UPF_sorted.bam"), path("*_UPR_sorted.bam"), path("*_SE_sorted.bam")

    script: 
    
    """
    dcsrsoft use ${params.softstack}
    module load ${params.samtools}
    new_name=${specimen_id}_${params.run_id}_${lane}

    samtools sort -@ ${task.cpus} -o \${new_name}_PE_sorted.bam -O bam ${PE}
    samtools sort -@ ${task.cpus} -o \${new_name}_UPF_sorted.bam -O bam ${UPF}
    samtools sort -@ ${task.cpus} -o \${new_name}_UPR_sorted.bam -O bam ${UPR}
    samtools sort -@ ${task.cpus} -o \${new_name}_SE_sorted.bam -O bam ${SE}
    """
}

/*
    merge sam files
*/

process merge_bams {
    label "merge_bams"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/2_mapped", mode: params.publish

    input: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path(PE), path(UPF), path(UPR), path(SE)
    output: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*_mapped.bam")

    script: 
    
    """
    dcsrsoft use ${params.softstack}
    module load ${params.samtools}

    new_name=${specimen_id}_${params.run_id}_${lane}

    samtools merge -f \
        -@ ${task.cpus} \
        -o \${new_name}_mapped.bam \
        ${PE} \
        ${UPF} \
        ${UPR} \
        ${SE}
    """

}

/*
    Mark duplicates & samtools index
    Needed to add indexing step to channel bam & bai together
*/

process markdup {
    label "markdup"
    tag "$specimen_id"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir/2_mapped", mode: params.publish, pattern: "*bam*"
    publishDir "$params.result_dir/QC/5_markedups", mode: params.publish, pattern: "*txt"

    input: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path(bam) 

    output: 
        tuple val(source_organism_id), val(specimen_id), val(lane), path("*_marked.bam"), path("*_marked.bam.bai"), emit: bams
        path("*_marked_stats.txt"), emit: stats
        path("*.tsv"), emit: tmptsv

    script:     
    """
    dcsrsoft use ${params.softstack}
    
    module load ${params.gcc} ${params.picard} ${params.samtools}
    
    new_name=${specimen_id}_${params.run_id}_${lane}

    picard MarkDuplicates -I ${bam} -O \${new_name}_marked.bam -M \${new_name}_marked_stats.txt

    samtools index \${new_name}_marked.bam

    echo -e "${source_organism_id}\t${specimen_id}\t${lane}\t\${new_name}_marked.bam" > \${new_name}_out_tmp.tsv
    """
}


/*
    Print output file
*/

process outtsv {
    label "outtsv"
    scratch params.scratch_dir // use $TMPDIR for I/O intensive tasks
    publishDir "$params.result_dir", mode: params.publish

     input: 
        path(tsvs)
    output: 
        path("*tsv")

    script:     
    """
    dcsrsoft use ${params.softstack}
    cat $tsvs > bamlist_${params.run_id}.tsv
    """
}


include { pretrim_qc_wf } from './modules/pretrimQC.nf'
include { posttrim_qc_wf } from './modules/posttrimQC.nf'
include { trim_stats } from './modules/trimQC.nf'
include { qualimap_wf} from './modules/qualimap.nf'
include { dups_multiQC} from './modules/dupsQC.nf'
include { sex_assign_wf} from './modules/sex_assignment.nf'

workflow {
    // load input files
    input = Channel
        .fromPath(params.input)
        .splitCsv( header: true, sep: ',' )        
        .map { row -> tuple( row.source_organism_id, row.specimen_id, row.lane, file(row.R1), file(row.R2) ) }
        .set { fastq_channel } 

    adapter_file = file(params.illumina_adapters)

    ref_index = file(params.ref_genome + ".{,amb,ann,bwt,pac,sa}")
    ref_fasta = file(params.ref_genome)

    lg_scaffolds = file(params.qualimap_bed)

    // process from raw reads to mapped & duplicate marked
    getRG(fastq_channel)
    merge_overlaping_reads(fastq_channel)
    trim_PE(fastq_channel, adapter_file)
    trim_SE(merge_overlaping_reads.out, adapter_file)
    bwa_input = trim_PE.out.trimmed_PE.join(trim_SE.out.trimmed_SE, by: [0, 1, 2])
    // bwa_input.view()
    bwa_map(bwa_input, ref_index)
    
    // join SAM & corresponding read group info
    rg_sam = bwa_map.out.aln.join(getRG.out, by: [0, 1, 2])
    read_groups(rg_sam)

    // sort SAM files
    sortbam(read_groups.out)
    
    // merge sorted PE, UPF, UPR and SE into single BAM
    merge_bams(sortbam.out)
    
    // mark duplicates and index
    markdup(merge_bams.out)
    
    // assign sex based on coverage
    sex_assign_wf(markdup.out.bams)

    // generate a tsv inventory of final BAM files
    bamlists = markdup.out.tmptsv.collect()
    bamlists.view()
    outtsv(bamlists)
    
    // QC
    pretrim_qc_wf(fastq_channel)
    posttrim_qc_wf(trim_PE.out.trimmed_PE, trim_SE.out.trimmed_SE)
    // trimstats_input = trim_PE.out.trimstats.collect().mix(trim_SE.out.trimstats.collect())
    trimstats_input = trim_PE.out.trimstats
    .mix(trim_SE.out.trimstats)
    .collect()
    // trimstats_input.view() // verify you get logs for paired, unpaired and assembled in 1 channel per library
    trim_stats(trimstats_input)
    qualimap_wf(markdup.out.bams, lg_scaffolds)
    dups_multiQC(markdup.out.stats.collect())
}