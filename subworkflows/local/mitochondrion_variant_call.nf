include { MERGE_STATS                             } from '../../modules/local/gatk/merge_mutect_stats.nf'
include { GET_CONTAMINATION                       } from '../../modules/local/haplocheckCLI/haplocheck_cli.nf'
include { LIFTOVER_VCF                            } from '../../modules/local/picard/liftover_vcf.nf'
include { MERGE_VCFS                              } from '../../modules/local/picard/merge_vcfs.nf'
include { FILTER_MUTECT_CALLS                     } from '../../modules/local/gatk/mitochondrial_variants_filter.nf'
include { LEFT_ALIGN_AND_TRIM_VARIANTS            } from '../../modules/local/gatk/left_align_and_trim_variants.nf'
include { SELECT_VARIANTS                         } from '../../modules/local/gatk/select_variants.nf'
include {
        CALL_VARIANTS as CALL_DEFAULT;
        CALL_VARIANTS as CALL_SHIFTED             } from '../../subworkflows/local/mutect2_variant_call.nf'



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

        LEFT_ALIGN_AND_TRIM_VARIANTS(
            params.mito_fasta,
            params.mito_index,
            params.mito_dict,
            FILTER_MUTECT_CALLS.out
        )

        SELECT_VARIANTS(LEFT_ALIGN_AND_TRIM_VARIANTS.out)

        GET_CONTAMINATION( SELECT_VARIANTS.out )

    emit:
        mutect_vcf = FILTER_MUTECT_CALLS.out                 // channel: [ val(sample_id), vcf, tbi
        contamination = GET_CONTAMINATION.out                // channel: [ val(sample_id), contam
        alignment = CALL_DEFAULT.out.alignment               // channel: [ val(sample_id), bam, bai
        alignment_metrics = CALL_DEFAULT.out.algn_metrics    // channel: [ val(sample_id), algn_metrics, theoretical_sensitivity
        alignment_wgs = CALL_DEFAULT.out.wgs_metrics         // channel: [ val(sample_id), wgs_metrics
        dup_metrics = CALL_DEFAULT.out.dup_metrics           // channel: [ val(sample_id), dup_metrics
}
