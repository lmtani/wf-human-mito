//
// Call mitochondrial variants for the given interval
//
include { BWA_ALIGN_FROM_UBAM as ALIGN_MITO  } from '../../modules/local/custom/bwa_align_from_ubam.nf'
include { GATK4_MUTECT2                      } from '../../modules/nf-core/gatk4/mutect2/main'
include { PICARD_COLLECTMULTIPLEMETRICS      } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_MARKDUPLICATES              } from '../../modules/nf-core/picard/markduplicates/main'
include { PICARD_SORTSAM                     } from '../../modules/nf-core/picard/sortsam/main'
include { SAMTOOLS_INDEX                     } from '../../modules/nf-core/samtools/index/main'


workflow CALL_VARIANTS {
    take:
        reads     // channel: [ val(meta), ubam ]
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
        ch_versions = Channel.empty()

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
        ch_versions = ch_versions.mix(ALIGN_MITO.out.versions)

        PICARD_MARKDUPLICATES(ALIGN_MITO.out.bam, fasta, fasta_fai)
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

        PICARD_SORTSAM(PICARD_MARKDUPLICATES.out.bam, "coordinate")
        ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions)

        SAMTOOLS_INDEX(PICARD_SORTSAM.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        PICARD_COLLECTMULTIPLEMETRICS(PICARD_SORTSAM.out.bam, fasta, fasta_fai)
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

        mutect_inputs = PICARD_SORTSAM.out.bam.join(SAMTOOLS_INDEX.out.bai).map{ it -> [ it[0], it[1], it[2], [] ] }
        GATK4_MUTECT2(mutect_inputs, fasta, fasta_fai, fasta_dict, [], [], [], [])
        ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    emit:
        alignment    = mutect_inputs                               // channel: [ val(meta), bam, bai ]
        dup_metrics  = PICARD_MARKDUPLICATES.out.metrics           // channel: [ val(meta), metrics ]
        mutect_stats = GATK4_MUTECT2.out.stats                     // channel: [ val(meta), stats ]
        pdf_metrics  = PICARD_COLLECTMULTIPLEMETRICS.out.pdf
        txt_metrics  = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        vcf          = GATK4_MUTECT2.out.vcf                       // channel: [ val(meta), vcf, tbi ]
        versions     = ch_versions                                 // channel: [ versions.yml ]
}
