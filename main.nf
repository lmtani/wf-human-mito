#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include {
    BWA_ALIGN_FROM_UBAM as ALIGN_SHIFTED_MITO;
    BWA_ALIGN_FROM_UBAM as ALIGN_STANDARD_MITO;
    BWA_ALIGN_FROM_UBAM;
    CREATE_ALL_SAMPLES_CSV;
    FILTER as INITIAL_FILTER;
    GET_CONTAMINATION;
    LIFTOVER_AND_COMBINE_VCFS;
    MERGE_STATS
    SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES} from './modules.nf'

include { FASTQ_TO_UBAM } from './modules/local/picard/fastq_to_sam.nf'

include { CREATE_JSON } from './modules/local/custom/create_sample_json.nf'

include { COLLECT_ALIGNMENT_METRICS } from './modules/local/picard/collect_alignment_summary_metrics.nf'

include { PRINT_READS } from './modules/local/gatk/print_reads.nf'
include { 
    CALL_MUTECT as CALL_MUTECT_SHIFTED;
    CALL_MUTECT as CALL_MUTECT_STANDARD } from './modules/local/gatk/mutect2.nf'

include { SELECT_MITO_READS } from './modules/local/picard/revert_sam.nf'
include { COLLECT_WGS_METRICS } from './modules/local/picard/collect_wgs_metrics.nf'

include { SORT_SAM; SORT_SAM as SORT_SHIFTED; SORT_SAM as SORT_DEFAULT } from './modules/local/picard/sort_sam.nf'

include { MARK_DUPLICATES as MARK_DUP_SHIFTED; MARK_DUPLICATES as MARK_DUP_DEFAULT } from './modules/local/picard/mark_duplicates.nf'


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
        ALIGN_STANDARD_MITO(
                reads,
                params.mito_fasta,
                params.mito_dict,
                params.mito_index,
                params.mito_amb,
                params.mito_ann,
                params.mito_bwt,
                params.mito_pac,
                params.mito_sa,
                params.mito_fake_alt
        )

        MARK_DUP_DEFAULT(ALIGN_STANDARD_MITO.out)

        SORT_DEFAULT(MARK_DUP_DEFAULT.out.bam)

        ALIGN_SHIFTED_MITO(
            reads,
            params.shifted_fasta,
            params.shifted_dict,
            params.shifted_index,
            params.shifted_amb,
            params.shifted_ann,
            params.shifted_bwt,
            params.shifted_pac,
            params.shifted_sa,
            params.mito_fake_alt
        )

        MARK_DUP_SHIFTED(ALIGN_SHIFTED_MITO.out)

        SORT_SHIFTED(MARK_DUP_SHIFTED.out.bam)

        COLLECT_ALIGNMENT_METRICS(
            SORT_DEFAULT.out,
            params.mito_fasta,
            params.mito_dict,
            params.mito_index
        )

        COLLECT_WGS_METRICS(
            SORT_DEFAULT.out,
            params.mito_fasta,
            300
        )

        CALL_MUTECT_STANDARD(
            SORT_DEFAULT.out,
            params.mito_fasta,
            params.mito_dict,
            params.mito_index,
            "standard",
            " -L chrM:576-16024 "
        )

        CALL_MUTECT_SHIFTED(
            SORT_SHIFTED.out,
            params.shifted_fasta,
            params.shifted_dict,
            params.shifted_index,
            "shifted",
            " -L chrM:8025-9144 "
        )

        ch2 = CALL_MUTECT_STANDARD.out.join(CALL_MUTECT_SHIFTED.out)
        LIFTOVER_AND_COMBINE_VCFS(
            ch2,
            params.mito_fasta,
            params.mito_index,
            params.mito_dict,
            params.shift_back_chain
        )

        MERGE_STATS(ch2)

        ch3 = LIFTOVER_AND_COMBINE_VCFS.out.join(MERGE_STATS.out)
        INITIAL_FILTER(
            ch3,
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
            INITIAL_FILTER.out
        )

        GET_CONTAMINATION(
            SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out
        )

        ch4 = INITIAL_FILTER.out.join(GET_CONTAMINATION.out)
        ch5 = ch4.join(ALIGN_STANDARD_MITO.out)
        ch6 = ch5.join(COLLECT_ALIGNMENT_METRICS.out)
        ch7 = ch6.join(COLLECT_WGS_METRICS.out)

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
