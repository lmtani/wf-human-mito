#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { separate_mitochondrion } from './subworkflows/local/separate_mitochondrion_reads.nf'
include { variant_call } from './subworkflows/local/mitochondrion_variant_call.nf'

def helpMessage(){
    log.info """
Human Mitochondrial Analysis Workflow'

Usage:
    nextflow run lmtani/wf-human-mito [options]

Script Options:
    --fastq        REGEX     Path to FASTQ directory. Quote is required. Ex: "/path/to/fastqs/*_R{1,2}*.fastq.gz" (required)
    --reference    FILE     Path to reference (GRCh38). BWA index need to be in same directory (required)
    --outdir      DIR     Path for output (default: $params.outdir)
"""
}


workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }

    if (!params.reference) {
        helpMessage()
        println("")
        println("`--reference` is required")
        exit 1
    }


    reads = Channel.fromFilePairs("${params.fastq}", glob: true)

    separate_mitochondrion(reads)

    variant_call(separate_mitochondrion.out)

    // CREATE_JSON(variant_call.out)
    // CREATE_ALL_SAMPLES_CSV(CREATE_JSON.out.collect())
}
