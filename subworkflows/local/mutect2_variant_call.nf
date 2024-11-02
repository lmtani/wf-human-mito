//
// Call mitochondrial variants for the given interval
//
include { BWA_ALIGN_FROM_UBAM as ALIGN_MITO  } from '../../modules/local/custom/bwa_align_from_ubam.nf'
include { GATK4_MUTECT2                      } from '../../modules/nf-core/gatk4/mutect2/main'
include { PICARD_COLLECTMULTIPLEMETRICS      } from '../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_MARKDUPLICATES              } from '../../modules/nf-core/picard/markduplicates/main'
include { PICARD_SORTSAM                     } from '../../modules/local/picard/sortsam/main'


workflow CALL_VARIANTS {
    take:
        reads             // channel: [ val(meta), ubam ]
        interval          // channel: val(interval)
        prefix            // channel: val(name)
        reference_genome  // channel: [ fasta, dict, index, amb, ann, bwt, pac, sa, alt, intervals ]

    main:
        // To gather all QC reports and versions for MultiQC
        ch_reports  = Channel.empty()
        ch_versions = Channel.empty()

        ALIGN_MITO(reads, reference_genome.dropRight(1))  // intervals are not needed for the alignment
        ch_versions = ch_versions.mix(ALIGN_MITO.out.versions)

        PICARD_MARKDUPLICATES(ALIGN_MITO.out.bam, reference_genome[0], reference_genome[2])
        ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

        PICARD_SORTSAM(PICARD_MARKDUPLICATES.out.bam, "coordinate")
        ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions)

        PICARD_COLLECTMULTIPLEMETRICS(PICARD_SORTSAM.out.bam.map{ it -> [ it[0], it[1] ] }, reference_genome[0], reference_genome[2])
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)

        mutect_inputs = PICARD_SORTSAM.out.bam.map{ it -> [ it[0], it[1], it[2], reference_genome[9] ] }
        GATK4_MUTECT2(mutect_inputs, reference_genome[0], reference_genome[2], reference_genome[1], [], [], [], [])
        ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions)

    emit:
        alignment      = mutect_inputs                               // channel: [ val(meta), bam, bai, intervals ]
        dup_metrics    = PICARD_MARKDUPLICATES.out.metrics           // channel: [ val(meta), metrics ]
        mutect_stats   = GATK4_MUTECT2.out.stats                     // channel: [ val(meta), stats ]
        pdf_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.pdf
        txt_metrics    = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
        vcf            = GATK4_MUTECT2.out.vcf                       // channel: [ val(meta), vcf, tbi ]
        versions       = ch_versions                                 // channel: [ versions.yml ]
        multiqc_rename = ALIGN_MITO.out.multiqc_rename               // channel: [ multiqc_rename.tsv ]
}
