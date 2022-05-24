#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include {
    GET_CONTAMINATION;
    SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES} from './modules.nf'


include { MERGE_STATS } from './modules/local/gatk/merge_mutect_stats.nf'
include { BWA_ALIGN_FROM_UBAM } from './modules.nf'
include { FASTQ_TO_UBAM } from './modules/local/picard/fastq_to_sam.nf'
include { LIFTOVER_VCF } from './modules/local/picard/liftover_vcf.nf'
include { MERGE_VCFS } from './modules/local/picard/merge_vcfs.nf'
include { PRINT_READS } from './modules/local/gatk/print_reads.nf'
include { FILTER_MUTECT_CALLS } from './modules/local/gatk/mitochondrial_variants_filter.nf'
include { SELECT_MITO_READS } from './modules/local/picard/revert_sam.nf'
include { SORT_SAM } from './modules/local/picard/sort_sam.nf'
include { CALL_VARIANTS as CALL_DEFAULT; CALL_VARIANTS as CALL_SHIFTED } from './subworkflows/local/mutect2_variant_call.nf'

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


workflow separate_mitochondrion {
    take: reads
    main:
        // Human Reference
        human_fasta = file("${params.reference}", type:'file', checkIfExists:true)
        human_dict = file("${human_fasta.getParent()}/${human_fasta.baseName}.dict", type:'file', checkIfExists:true)
        human_index = file("${params.reference}.fai", type:'file', checkIfExists:true)
        human_amb = file("${params.reference}.64.amb", type:'file', checkIfExists:true)
        human_ann = file("${params.reference}.64.ann", type:'file', checkIfExists:true)
        human_bwt = file("${params.reference}.64.bwt", type:'file', checkIfExists:true)
        human_sa = file("${params.reference}.64.sa", type:'file', checkIfExists:true)
        human_pac = file("${params.reference}.64.pac", type:'file', checkIfExists:true)
        human_alt = file("${params.reference}.64.alt", type:'file', checkIfExists:true)

        FASTQ_TO_UBAM(reads)
        BWA_ALIGN_FROM_UBAM(
            FASTQ_TO_UBAM.out,
            human_fasta,
            human_dict,
            human_index,
            human_amb,
            human_ann,
            human_bwt,
            human_pac,
            human_sa,
            human_alt
        )

        SORT_SAM(BWA_ALIGN_FROM_UBAM.out)

        PRINT_READS(
            SORT_SAM.out,
            human_fasta,
            human_index,
            human_dict,
        )

        SELECT_MITO_READS(
            PRINT_READS.out,
            human_fasta,
            human_dict,
            human_index,
        )
    emit:
        SELECT_MITO_READS.out
}


workflow variant_call {
    take:
        reads
    main:
        CALL_DEFAULT(reads, " -L chrM:576-16024 ", "standard")
        CALL_SHIFTED(reads, " -L chrM:8025-9144 ", "shifted")

        LIFTOVER_VCF(
            CALL_SHIFTED.out.variants,
            params.mito_fasta,
            params.mito_index,
            params.mito_dict,
            params.shift_back_chain,
        )

        MERGE_VCFS(CALL_DEFAULT.out.variants.join(LIFTOVER_VCF.out))
        MERGE_STATS(CALL_DEFAULT.out.mutect_stats.join(CALL_SHIFTED.out.mutect_stats))

        FILTER_MUTECT_CALLS(
            MERGE_VCFS.out.join(MERGE_STATS.out),
            params.mito_fasta,
            params.mito_index,
            params.mito_dict,
            params.blacklist,
            params.blacklist_index
        )

        SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES(
            params.mito_fasta,
            params.mito_index,
            params.mito_dict,
            FILTER_MUTECT_CALLS.out
        )

        GET_CONTAMINATION(
            SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out
        )

        ch4 = FILTER_MUTECT_CALLS.out.join(GET_CONTAMINATION.out)
        ch5 = ch4.join(CALL_DEFAULT.out.alignment)
        ch6 = ch5.join(CALL_DEFAULT.out.algn_metrics)
        ch7 = ch6.join(CALL_DEFAULT.out.wgs_metrics)

    emit:
        ch7
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
    reads.view()

    separate_mitochondrion(reads)

    variant_call(separate_mitochondrion.out)

    // CREATE_JSON(variant_call.out)
    // CREATE_ALL_SAMPLES_CSV(CREATE_JSON.out.collect())
}
