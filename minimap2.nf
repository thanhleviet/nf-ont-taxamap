#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process DEHUMAN {
    publishDir "${params.outdir}", mode: "copy"
    
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


process MINIMAP2 {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {sample_id}
    
    cpus 32

    memory '64.GB'
    
    maxForks 2
    
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.paf.gz")

    script:
    """
    minimap2 -cx map-ont -t $task.cpus $params.ref $reads | pigz - > ${sample_id}.paf.gz
    """
}


process CLEAN_PAF {
    publishDir "${params.outdir}", mode: "copy"
    
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
    paf2taxid.py -o ${sample_id} --score ${params.score} ${paf}
    """
}

ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
                  .map {it -> tuple(it.simpleName, it)}

log.info "Ref database: " + params.ref
log.info "Outdir: " + params.outdir
log.info "Score: " + params.score
params.score = null

workflow {
    DEHUMAN(ch_reads)
    MINIMAP2(DEHUMAN.out)
    PARSING_PAF(MINIMAP2.out)
    CLEAN_PAF(MINIMAP2.out)
}