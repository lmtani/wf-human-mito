//
// Call mitochondrial variants, using shifted approach
// described in https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/#mtdna-calling-pipeline-for-single-samples
//
include { MERGE_STATS                             } from '../../modules/local/gatk/merge_mutect_stats.nf'
include { HAPLOCHECK                              } from '../../modules/local/haplocheck/haplocheck.nf'
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
        reads  // channel: [ val(meta), ubam ]
    main:
        CALL_DEFAULT(
            reads,
            " -L chrM:576-16024 ",
            "standard",
            params.genome.mito_fasta,
            params.genome.mito_dict,
            params.genome.mito_index,
            params.genome.mito_amb,
            params.genome.mito_ann,
            params.genome.mito_bwt,
            params.genome.mito_pac,
            params.genome.mito_sa,
            params.genome.mito_fake_alt,
        )
        CALL_SHIFTED(
            reads,
            " -L chrM:8025-9144 ",
            "shifted",
            params.genome.shifted_fasta,
            params.genome.shifted_dict,
            params.genome.shifted_index,
            params.genome.shifted_amb,
            params.genome.shifted_ann,
            params.genome.shifted_bwt,
            params.genome.shifted_pac,
            params.genome.shifted_sa,
            params.genome.mito_fake_alt,
        )

    //     LIFTOVER_VCF(
    //         CALL_SHIFTED.out.variants,
    //         params.genome.mito_fasta,
    //         params.genome.mito_index,
    //         params.genome.mito_dict,
    //         params.genome.shift_back_chain,
    //     )

    //     MERGE_VCFS(CALL_DEFAULT.out.variants.join(LIFTOVER_VCF.out))
    //     MERGE_STATS(CALL_DEFAULT.out.mutect_stats.join(CALL_SHIFTED.out.mutect_stats))

    //     FILTER_MUTECT_CALLS(
    //         MERGE_VCFS.out.join(MERGE_STATS.out),
    //         params.genome.mito_fasta,
    //         params.genome.mito_index,
    //         params.genome.mito_dict,
    //         params.genome.blacklist,
    //         params.genome.blacklist_index
    //     )

    //     LEFT_ALIGN_AND_TRIM_VARIANTS(
    //         params.genome.mito_fasta,
    //         params.genome.mito_index,
    //         params.genome.mito_dict,
    //         FILTER_MUTECT_CALLS.out
    //     )

    //     SELECT_VARIANTS(LEFT_ALIGN_AND_TRIM_VARIANTS.out)

    //     HAPLOCHECK(SELECT_VARIANTS.out)

    // emit:
    //     mutect_vcf        = FILTER_MUTECT_CALLS.out        // channel: [ val(sample_id), vcf, tbi ]
    //     contamination     = HAPLOCHECK.out.   txt          // channel: [ val(sample_id), contam ]
        alignment         = CALL_DEFAULT.out.alignment     // channel: [ val(sample_id), bam, bai ]
        // alignment_metrics = CALL_DEFAULT.out.algn_metrics  // channel: [ val(sample_id), algn_metrics, theoretical_sensitivity ]
        // alignment_wgs     = CALL_DEFAULT.out.wgs_metrics   // channel: [ val(sample_id), wgs_metrics ]
        dup_metrics       = CALL_DEFAULT.out.dup_metrics   // channel: [ val(sample_id), dup_metrics ]
}
