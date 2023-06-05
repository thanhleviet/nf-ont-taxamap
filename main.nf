#!/usr/bin/env nextflow
include { paramsHelp } from 'plugin/nf-validation'
include { paramsSummaryLog } from 'plugin/nf-validation'
include { scanFolder } from './lib/utils.nf'

log.info paramsSummaryLog(workflow)

if (params.help) {
    log.info paramsHelp("nextflow run thanhleviet/ont-meta-minimap2 --input /path/to/a/folder/of/barcodes")
    exit 0
}

process CONCATENATE {
    
    tag {sample_id}
    
    conda './env/conda-env.yml'
    // conda 'epi2melabs::fastcat'

    cpus 4

    input:
        tuple val(sample_id), path(barcode_path)

    output:
        tuple val(sample_id), path("${sample_id}.fastq.gz")

    script:
    """
    fastcat ${barcode_path} | pigz -p ${task.cpus} - >  ${sample_id}.fastq.gz
    """
    
    stub:
    """
    fastcat --help
    touch ${sample_id}.fastq.gz
    """
}

process DEHUMAN {
    
    tag {sample_id}
    
    conda './env/conda-env.yml'

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
    stub:
    """
    minimap2 --version
    touch ${sample_id}.non_host.fastq.gz
    """
}

process FILTER_BY_LENGTH {
    
    tag {sample_id}
    conda './env/conda-env.yml'

    cpus 4

    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path("${sample_id}.filtered.${params.min_read_length}.fastq.gz"), emit: reads
        tuple val(sample_id), path("*.json"), emit: stats

    script:
    """
    seqkit -j $task.cpus seq -m ${params.min_read_length} ${reads} | pigz - > ${sample_id}.filtered.${params.min_read_length}.fastq.gz
    nanoq -s -vvv -j -r ${sample_id}.stats.pre_filtered.json -i ${reads}
    nanoq -s -vvv -j -r ${sample_id}.stats.post_filtered.json -i ${sample_id}.filtered.${params.min_read_length}.fastq.gz
    """
    stub:
    """
    seqkit -h
    touch ${sample_id}.filtered.${params.min_read_length}.fastq.gz
    touch ${sample_id}.stats.pre_filtered.json
    touch ${sample_id}.stats.post_filtered.json 
    """
}

process DEHUMAN2 {
    
    tag {sample_id}
    conda './env/conda-env.yml'

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
    stub:
    """
    minimap2 --version
    samtools --version
    touch ${sample_id}.non_host.fastq.gz
    touch ${sample_id}.flagstat
    """
}

process MINIMAP2_PAF {
    
    tag {sample_id}
    
    conda './env/conda-env.yml'

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
    stub:
    """
    minimap2 -h
    touch ${sample_id}.paf.gz
    """

}

process AMR {
    
    tag {sample_id}
    
    conda './env/conda-env.yml'

    cpus 32

    memory '64.GB'
    
    // maxForks 2
    
    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.*")

    script:
    """
    kma -i ${reads} -cge -ont -bcNano -md 1.0 -t ${task.cpus} -t_db ${params.amr_db} -o ${sample_id}
    """
    stub:
    """
    kma -v
    touch ${sample_id}.{res,sam,map}
    """
}


process MINIMAP2_SAM {
    
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
    stub:
    """
    minimap2 -v
    touch ${sample_id}.bam
    """
}

process CLEAN_PAF {
    
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
    stub:
    """
    echo "Checked!"
    touch ${sample_id}.cleaned.paf.gz
    """
}

process SAM2LCA {
    
    tag {sample_id}
    
    cpus 12

    input:
    tuple val(sample_id), path(bam)
    output:
    path("*.{csv,json,bam,report}")

    script:
    """
    samtools index ${bam}
    sam2lca -d ${params.sam2lca_db} analyze -b -p ${task.cpus} ${bam}
    """
    stub:
    """
    samtools -v
    """
}

process PARSING_PAF {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {sample_id}
    
    conda './env/conda-report.yml'
    
    cpus 4

    input:
        tuple val(sample_id), path(paf)
    output:
        path("*.csv")
    script:
    """
    paf2taxid.py -o ${sample_id} --score ${params.score} --min_read_length ${params.min_read_length} ${paf}
    """
    stub:
    """
    echo "Checked!"
    touch ${sample_id}.csv
    """
}

// ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
//                   .map {it -> tuple(it.simpleName, it)}

// log.info "Ref database: " + params.ref
// log.info "Ref_Sam database: " + params.ref_sam
// log.info "Outdir: " + params.outdir
// log.info "Score: " + params.score


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

def folderContents = scanFolder(params.input)
def input_type = folderContents[0][0]
def pattern = folderContents[0][1]
ch_input = Channel.fromPath("${params.input}/${pattern}", type: input_type)
                    .map {it -> tuple(it.simpleName, it)}


workflow {

    if (input_type == "dir") {
        ch_reads = CONCATENATE(ch_input)
    } else {
        ch_reads = ch_input
    }
    
    FILTER_BY_LENGTH(ch_reads)
    
    if (!params.skip_dehuman) {
        DEHUMAN2(FILTER_BY_LENGTH.out.reads)
        ch_dehuman_reads = DEHUMAN2.out.fastq
    } else {
        ch_dehuman_reads = FILTER_BY_LENGTH.out.reads
    }

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