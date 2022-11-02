#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { separate_mitochondrion      } from './subworkflows/local/separate_mitochondrion_reads.nf'
include { variant_call                } from './subworkflows/local/mitochondrion_variant_call.nf'
include { make_report                 } from './subworkflows/local/make_report.nf'
include { MULTIQC                     } from './modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/local/custom/dumpsoftwareversions/main'

def helpMessage(){
    log.info """
    Human Mitochondrial Analysis Workflow'

    Usage:
        nextflow run lmtani/wf-human-mito [options]

    Script Options:
        --fastq             REGEX     Path to FASTQ directory. Quote is required. Ex: "/path/to/fastqs/*_R{1,2}*.fastq.gz" (optional if alignments provided)
        --alignments        REGEX     Path to the directory with alignments in BAM or CRAM format. (optional if fastq provided)
        --reference         FILE      Path to reference (GRCh38). BWA index need to be in same directory (required)
        --outdir            DIR       Path for output (default: $params.outdir)
        --restore_hardclips BOOLEAN   When true, restores reads and qualities of records with hard-clips containing XB and XQ tags (useful when your inputs are alignments)
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
    ch_versions = Channel.empty()
    ch_multiqc_files =  Channel.empty()

    if (params.fastq) {
        reads = Channel.fromFilePairs("${params.fastq}", glob: true)
    }
    if (params.alignments) {
        alignments = Channel.fromFilePairs("${params.alignments}", glob: true, flat: true).map( it ->
            [ [id: it[0]], it[2], it[1] ]
        )
    }

    separate_mitochondrion(reads, alignments, params.restore_hardclips)
    ch_versions = ch_versions.mix(separate_mitochondrion.out.versions)

    variant_call(separate_mitochondrion.out.bam)
    ch_versions = ch_versions.mix(variant_call.out.versions)
    ch_reports = ch_versions.mix(variant_call.out.reports)

    make_report(variant_call.out.contamination)
    
    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
    ch_multiqc_files = ch_multiqc_files.mix(
        variant_call.out.multiqc_rename.unique().collectFile(name: 'rename_samples.tsv').collect(), 
        ch_reports.collect().ifEmpty([]), 
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    )
    
    MULTIQC(ch_multiqc_files.collect(), [], [], [])
}
