//
// Call mitochondrial variants
//


include { BWA_ALIGN_FROM_UBAM           } from '../../modules.nf'
include { PRINT_READS                   } from '../../modules/local/gatk/print_reads.nf'
include { CALL_MUTECT                   } from '../../modules/local/gatk/mutect2.nf'
include { COLLECT_WGS_METRICS           } from '../../modules/local/picard/collect_wgs_metrics.nf'
include { COLLECT_ALIGNMENT_METRICS     } from '../../modules/local/picard/collect_alignment_summary_metrics.nf'
include { MARK_DUPLICATES               } from '../../modules/local/picard/mark_duplicates.nf'
include { SORT_SAM                      } from '../../modules/local/picard/sort_sam.nf'


workflow CALL_VARIANTS {
    take:
        reads     // channel: [ val(sample_id), reads_input  ]
        interval  // channel: [ val(interval) ]
        prefix    // channel: [ val(name) ]

    main:
        BWA_ALIGN_FROM_UBAM(
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

        MARK_DUPLICATES(BWA_ALIGN_FROM_UBAM.out)

        SORT_SAM(MARK_DUPLICATES.out.bam)

        COLLECT_ALIGNMENT_METRICS(
            SORT_SAM.out,
            params.mito_fasta,
            params.mito_dict,
            params.mito_index
        )

        COLLECT_WGS_METRICS(
            SORT_SAM.out,
            params.mito_fasta,
            300
        )

        CALL_MUTECT(
            SORT_SAM.out,
            params.mito_fasta,
            params.mito_dict,
            params.mito_index,
            prefix,
            interval
        )

    emit:
        variants             = CALL_MUTECT.out.vcf            // channel: [ val(sample_id), vcf, tbi ]
        mutect_stats         = CALL_MUTECT.out.stats          // channel: [ val(sample_id), stats ]
        wgs_metrics          = COLLECT_WGS_METRICS.out        // channel: [ val(sample_id), theoretical_sensibility, metrics ]
        algn_metrics         = COLLECT_ALIGNMENT_METRICS.out  // channel: [ val(sample_id), metrics ]
        alignment            = SORT_SAM.out                   // channel: [ val(sample_id), bam, bai ]
}