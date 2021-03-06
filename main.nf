#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { separate_mitochondrion } from './subworkflows/local/separate_mitochondrion_reads.nf'
include { variant_call           } from './subworkflows/local/mitochondrion_variant_call.nf'
include { make_report            } from './subworkflows/local/make_report.nf'


def helpMessage(){
    log.info """
    Human Mitochondrial Analysis Workflow'

    Usage:
        nextflow run lmtani/wf-human-mito [options]

    Script Options:
        --fastq        REGEX     Path to FASTQ directory. Quote is required. Ex: "/path/to/fastqs/*_R{1,2}*.fastq.gz" (optional if alignments provided)
        --alignments   REGEX     Path to the directory with alignments in BAM or CRAM format. (optional if fastq provided)
        --reference    FILE      Path to reference (GRCh38). BWA index need to be in same directory (required)
        --outdir       DIR       Path for output (default: $params.outdir)
    """.stripIndent()
}


workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if ((!params.fastq) && (!params.alignments)){
        helpMessage()
        println("")
        println("`--fastq` or --alignments required")
        exit 1
    }

    if (!params.reference) {
        helpMessage()
        println("")
        println("`--reference` is required")
        exit 1
    }

    reads = Channel.empty()
    alignments = Channel.empty()

    if (params.fastq) {
        reads = Channel.fromFilePairs("${params.fastq}", glob: true)
    }
    if (params.alignments) {
        alignments = Channel.fromFilePairs("${params.alignments}", glob: true, flat: true)
    }

    separate_mitochondrion(reads, alignments)

    variant_call(separate_mitochondrion.out)

    make_report(
        variant_call.out.contamination,
        variant_call.out.alignment_metrics,
        variant_call.out.alignment_wgs,
        variant_call.out.dup_metrics
    )
}
