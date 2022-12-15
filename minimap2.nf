#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DEHUMAN {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 32

    memory '64.GB'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.non_host.fastq.gz")

    script:
    """
    minimap2 -ax map-ont -t $task.cpus $params.human_ref $reads | samtools fastq -f 4 | pigz - > ${sample_id}.non_host.fastq.gz
    """
}

process DEHUMAN2 {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 32

    memory '64.GB'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.non_host.fastq.gz"), emit: fastq
        path ("${sample_id}.flagstat"), emit: flagstat

    script:
    """
    minimap2 -ax map-ont -t $task.cpus $params.human_ref $reads | samtools view -Sb - > mapped.bam
    samtools flagstat mapped.bam > ${sample_id}.flagstat
    samtools fastq -f 4 mapped.bam | pigz - > ${sample_id}.non_host.fastq.gz
    """
}

process MINIMAP2_PAF {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 32

    memory '64.GB'
    
    // maxForks 2
    
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.paf.gz")

    script:
    """
    minimap2 -cx map-ont -t $task.cpus $params.ref $reads --secondary=no | pigz - > ${sample_id}.paf.gz
    """
}


process AMR {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 32

    memory '64.GB'
    
    // maxForks 2
    
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.*")

    script:
    """
    kma -i ${reads} -ont -bcNano -md 1.0 -t ${task.cpus} -1t1 -t_db ${params.amr_db} -o ${sample_id}
    """
}


process MINIMAP2_SAM {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 32

    memory '64.GB'
    
    // maxForks 2
    
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    minimap2 -ax map-ont --split-prefix=tmp -N ${params.secondary_reads} -t $task.cpus $params.ref_sam $reads | samtools view -Sb - | samtools sort -@ $task.cpus -o ${sample_id}.bam -
    """
}

process CLEAN_PAF {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 4

    input:
        tuple val(sample_id), path(paf)
    output:
        path("${sample_id}.cleaned.paf.gz")
    
    shell:
    '''
        zcat !{sample_id}.paf.gz | awk 'OFS="," {print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$15,$17}' | sed 's/kraken:taxid|//g' | csvtk -H mutate2 -n cov -e '($4-$3)/$2*100' | csvtk -H mutate2 -n iden -e '$9/$10*100' | csvtk -H mutate2 -n score -e '2*(($14*$15)/($14+$15)*100)'| pigz - >  !{sample_id}.cleaned.paf.gz
    '''
}

process SAM2LCA {
    publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", mode: "copy"
    
    tag {sample_id}
    
    cpus 12

    input:
    tuple val(sample_id), path(bam)
    output:
    path("*.{csv,json,bam,report}")

    script:
    """
    hostname
    samtools index ${bam}
    sam2lca -d ${params.sam2lca_db} analyze -b -p ${task.cpus} ${bam}
    """
}

process PARSING_PAF {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {sample_id}
    
    cpus 4

    input:
        tuple val(sample_id), path(paf)
    output:
        path("*.csv")
    script:
    """
    paf2taxid.py -o ${sample_id} --score ${params.score} --min_read_length ${params.min_read_length} ${paf}
    """
}

ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
                  .map {it -> tuple(it.simpleName, it)}

log.info "Ref database: " + params.ref
log.info "Ref_Sam database: " + params.ref_sam
log.info "Outdir: " + params.outdir
log.info "Score: " + params.score

workflow wf_MINIMAP2_PAF {
    take:
        ch_reads
    main:
        MINIMAP2_PAF(ch_reads)
        PARSING_PAF(MINIMAP2_PAF.out)
        CLEAN_PAF(MINIMAP2_PAF.out)
}

workflow wf_MINIMAP2_SAM {
    take:
        ch_reads
    main:
        MINIMAP2_SAM(ch_reads)
        SAM2LCA(MINIMAP2_SAM.out)
}

workflow {
    DEHUMAN2(ch_reads)
    ch_dehuman_reads = DEHUMAN2.out.fastq
    if (!params.skip_amr) {
        AMR(ch_dehuman_reads)
    }
    if (!params.skip_paf) {
        wf_MINIMAP2_PAF(ch_dehuman_reads)
    }
    if (!params.skip_sam) {
        wf_MINIMAP2_SAM(ch_dehuman_reads)
    }
}