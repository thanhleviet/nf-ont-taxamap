#!/usr/bin/env nextflow

nextflow.enable.dsl=2


params.list = "./blood_infectious_species.txt"

ch_list = Channel.fromPath(params.list)
                 .splitCsv(header: false)
                 .map {it -> it[0]}

process FILTER {
    publishDir  params.outdir, mode: "copy"
    
    tag {species}
    
    cpus 2

    input:
        val(species)
    output:

    script:
        """
        seqkit grep -n -r -p '${species}' ${params.fasta} > '${species}.fna'
        """
}

workflow {
    FILTER(ch_list)
}