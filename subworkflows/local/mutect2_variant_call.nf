//
// Call mitochondrial variants for the given interval
//
include { BWA_ALIGN_FROM_UBAM as ALIGN_MITO  } from '../../modules/local/custom/bwa_align_from_ubam.nf'
include { CALL_MUTECT                        } from '../../modules/local/gatk/mutect2.nf'
include { COLLECT_WGS_METRICS                } from '../../modules/local/picard/collect_wgs_metrics.nf'
include { COLLECT_ALIGNMENT_METRICS          } from '../../modules/local/picard/collect_alignment_summary_metrics.nf'
include { MARK_DUPLICATES                    } from '../../modules/local/picard/mark_duplicates.nf'
include { SORT_SAM                           } from '../../modules/local/picard/sort_sam.nf'  // TODO: use PICARD_SORTSAM instead


workflow CALL_VARIANTS {
    take:
        reads     // channel: [ val(sample_id), ubam ]
        interval  // channel: val(interval)
        prefix    // channel: val(name)
        fasta
        fasta_dict
        fasta_fai
        fasta_amb
        fasta_ann
        fasta_bwt
        fasta_pac
        fasta_sa
        fasta_alt

    main:
        ALIGN_MITO(
                reads,
                fasta,
                fasta_dict,
                fasta_fai,
                fasta_amb,
                fasta_ann,
                fasta_bwt,
                fasta_pac,
                fasta_sa,
                fasta_alt
        )

        MARK_DUPLICATES(ALIGN_MITO.out)

        SORT_SAM(MARK_DUPLICATES.out.bam)

        COLLECT_ALIGNMENT_METRICS(SORT_SAM.out, fasta, fasta_dict, fasta_fai)

        COLLECT_WGS_METRICS(SORT_SAM.out, fasta, 300)  //TODO: parse READ_LENGTH value

        CALL_MUTECT(SORT_SAM.out, fasta, fasta_dict, fasta_fai, prefix, interval)

    emit:
        variants             = CALL_MUTECT.out.vcf            // channel: [ val(sample_id), vcf, tbi ]
        mutect_stats         = CALL_MUTECT.out.stats          // channel: [ val(sample_id), stats ]
        wgs_metrics          = COLLECT_WGS_METRICS.out        // channel: [ val(sample_id), theoretical_sensibility, metrics ]
        algn_metrics         = COLLECT_ALIGNMENT_METRICS.out  // channel: [ val(sample_id), metrics ]
        alignment            = SORT_SAM.out                   // channel: [ val(sample_id), bam, bai ]
        dup_metrics          = MARK_DUPLICATES.out.metrics    // channel: [ val(sample_id), metrics ]
}
