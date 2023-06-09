#!/usr/bin/env nextflow
import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.attribute.BasicFileAttributes

include { paramsHelp } from 'plugin/nf-validation'
include { paramsSummaryLog } from 'plugin/nf-validation'
include { scanFolder; getMostRecentXlsxFile } from './lib/utils.nf'

log.info paramsSummaryLog(workflow)

if (params.help) {
    log.info paramsHelp("nextflow run thanhleviet/ont-meta-minimap2 --input /path/to/a/folder/of/barcodes")
    exit 0
}

process CONCATENATE {

    tag {sample_id}

    conda './env/conda-env.yml'

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

    cpus 6

    memory '12.GB'

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.non_host.fastq.gz"), emit: fastq
        path ("${sample_id}.flagstat"), emit: flagstat

    script:
    args = task.ext.args ?: ''
    """
    minimap2 -ax map-ont $args -t $task.cpus $params.human_ref $reads | samtools view -Sb - > mapped.bam
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

    cpus 14

    memory '64.GB'


    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}.paf.gz")

    script:
    args = task.ext.args ?: ''

    """
    minimap2 -cx map-ont $args -t $task.cpus $params.ref $reads --secondary=no | pigz - > ${sample_id}.paf.gz
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
    publishDir "${params.outdir}/report", mode: "copy"

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

process CREATE_EXCEL_REPORT {
    publishDir "${params.outdir}", mode: "copy"
    conda './env/conda-r.yml'

    cpus 4

    input:
        path("*")

    output:
        path("*.xlsx")

    script:
    // Creates a new `java.util.Date` object and formats it as a string with the pattern `yyyy-MM-dd_HH`
    _this_run = new java.util.Date().format('yyyy-MM-dd_HH')
    // Sets the variable `run_name` to either the value of the `params.run_name` parameter,
    // or to the formatted date string if `params.run_name` is null or undefined
    run_name = params.run_name ?: workflow.runName
    """
        Rscript --vanilla ${projectDir}/bin/create_excel_report.R \$PWD "report_${run_name}"
    """
    stub:
    """
    touch report.xlsx
    """
}


workflow wf_MINIMAP2_PAF {
    take:
        ch_reads
    main:
        MINIMAP2_PAF(ch_reads)
        PARSING_PAF(MINIMAP2_PAF.out)
        CLEAN_PAF(MINIMAP2_PAF.out)
    emit:
    // Maps the output of the process `PARSING_PAF.out` to a list of files that match the pattern `*_metaphlan_report.csv`
    report = PARSING_PAF.out.map { it ->
        def metaphlan_report = it.grep(~/.*_metaphlan_report\.csv$/)
        if (metaphlan_report) {
            metaphlan_report[0]
        } else {
            null
        }
    }
}

workflow wf_MINIMAP2_SAM {
    take:
        ch_reads
    main:
        MINIMAP2_SAM(ch_reads)
        SAM2LCA(MINIMAP2_SAM.out)
}

// Scans the folder specified in the `params.input` parameter and retrieves the type and pattern of the files in the folder
def folderContents = scanFolder(params.input)
def input_type = folderContents[0][0]
def pattern = folderContents[0][1]

// Creates a channel `ch_input` from the files in the folder that match the pattern,
// and maps each file to a tuple containing its simple name and its full path
ch_input = Channel.fromPath("${params.input}/${pattern}", type: input_type)
                .map {it -> tuple(it.simpleName, it)}


workflow {

    ch_report = Channel.empty()

    if (input_type == "dir") {
        ch_reads = CONCATENATE(ch_input)
    } else {
        ch_reads = ch_input
    }

    FILTER_BY_LENGTH(ch_reads)

    ch_report = ch_report.mix(FILTER_BY_LENGTH.out.stats.map {it -> it[1]}).ifEmpty([])

    if (!params.skip_dehuman) {
        DEHUMAN2(FILTER_BY_LENGTH.out.reads)
        ch_dehuman_reads = DEHUMAN2.out.fastq
        ch_report = ch_report.mix(DEHUMAN2.out.flagstat).ifEmpty([])
    } else {
        ch_dehuman_reads = FILTER_BY_LENGTH.out.reads
    }

    if (!params.skip_amr) {
        AMR(ch_dehuman_reads)
        // AMR.out.map{ it -> it[1].grep(~/.*\.res$/)}.view()

        // Mixes the channel `ch_report` with the output of the process `AMR.out`
        // and maps the output to a list of `.res` files
        ch_report = ch_report.mix(
            AMR.out.map {it ->
                def resFiles = it[1].grep(~/.*\.res$/)
                if (resFiles) {
                    resFiles[0]
                } else {
                    null
                }
            }
        ).ifEmpty([])
    }

    if (!params.skip_paf) {
        wf_MINIMAP2_PAF(ch_dehuman_reads)
    }

    // Mixes the channel `ch_report` with the output of the process `wf_MINIMAP2_PAF.out.report`
    // and if the resulting channel is empty, sets it to an empty list
    ch_report = ch_report.mix(wf_MINIMAP2_PAF.out.report).ifEmpty([])
    // Not used for now
    // if (!params.skip_sam) {
    //     wf_MINIMAP2_SAM(ch_dehuman_reads)
    // }

    // Collects all items from the channel `ch_report` and returns them as a list
    ch_report = ch_report.collect()
    CREATE_EXCEL_REPORT(ch_report)
}

workflow.onComplete {
    report = getMostRecentXlsxFile(params.outdir.toString())
    monitor_email = "thanh.le-viet@quadram.ac.uk"
    if (params.email) {
        email = params.email + "," + monitor_email
    } else {
        email = monitor_email
    }

    if (workflow.success) {
        def msg = """\
        Pipeline execution summary
        ---------------------------
        Run Name    : ${workflow.runName}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail {
            to email
            subject "BART pipeline completed!"
            body msg
            attach report
            }
    } else {
        nextflow_log = "${workflow.launchDir}/.nextflow.log"
    sendMail {
            to monitor_email
            subject "BART pipeline failed!" + workflow.runName
            body "Please check the log file for more details."
            attach nextflow_log
            }
    }
}
