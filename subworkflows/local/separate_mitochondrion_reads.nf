//
// Aligns paired FASTQ reads into the human reference and selects
// reads from mitochondrial genome.
//
include { BWA_ALIGN_FROM_UBAM as ALIGN_RAW_READS  } from '../../modules/local/custom/bwa_align_from_ubam.nf'
include { GATK4_FASTQTOSAM                        } from '../../modules/nf-core/gatk4/fastqtosam/main'
include { PICARD_SORTSAM                          } from '../../modules/nf-core/picard/sortsam/main'
include { PRINT_READS                             } from '../../modules/local/gatk/print_reads.nf'
include { SAMTOOLS_INDEX                          } from '../../modules/nf-core/samtools/index/main'
include { SELECT_MITO_READS                       } from '../../modules/local/picard/revert_sam.nf'


workflow separate_mitochondrion {
    take:
        reads             // channel: [ val(meta), [ path(fq_r1), path(fq_r2) ] ]
        alignments        // channel: [ val(meta), [ path(bam), path(bai) ] ]
        restore_hardclips // channel: [ val(boolean) ]
    main:
        // Human Reference
        fasta = file("${params.reference}", type:'file', checkIfExists:true)
        dict  = file("${fasta.getParent()}/${fasta.baseName}.dict", type:'file', checkIfExists:true)
        index = file("${params.reference}.fai", type:'file', checkIfExists:true)
        amb   = file("${params.reference}.64.amb", type:'file', checkIfExists:true)
        ann   = file("${params.reference}.64.ann", type:'file', checkIfExists:true)
        bwt   = file("${params.reference}.64.bwt", type:'file', checkIfExists:true)
        sa    = file("${params.reference}.64.sa", type:'file', checkIfExists:true)
        pac   = file("${params.reference}.64.pac", type:'file', checkIfExists:true)
        alt   = file("${params.reference}.64.alt", type:'file', checkIfExists:true)
        human_reference_genome = [ fasta, dict, index, amb, ann, bwt, sa, pac, alt]

        ch_versions = Channel.empty()

        // Create channel with 'meta' info
        sample = reads.map( it -> { [ [id:it[0], single_end:false], [ file(it[1][0]), file(it[1][1]) ] ] } )

        GATK4_FASTQTOSAM(sample)
        ch_versions = ch_versions.mix(GATK4_FASTQTOSAM.out.versions)

        ALIGN_RAW_READS(GATK4_FASTQTOSAM.out.bam, human_reference_genome)
        ch_versions = ch_versions.mix(ALIGN_RAW_READS.out.versions)

        PICARD_SORTSAM(ALIGN_RAW_READS.out.bam, "coordinate")
        ch_versions = ch_versions.mix(PICARD_SORTSAM.out.versions)

        SAMTOOLS_INDEX(PICARD_SORTSAM.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        // Join the mapped bam + bai paths by their keys
        bam_sorted_indexed = PICARD_SORTSAM.out.bam.join(SAMTOOLS_INDEX.out.bai)

        // Add provided alignments to the reads channel
        mito_reads_ch = bam_sorted_indexed.mix(alignments)

        PRINT_READS(mito_reads_ch, fasta, index, dict)
        ch_versions = ch_versions.mix(PRINT_READS.out.versions)

        SELECT_MITO_READS(PRINT_READS.out.bam, fasta, dict, index, restore_hardclips)
        ch_versions = ch_versions.mix(SELECT_MITO_READS.out.versions)

    emit:
        bam      = SELECT_MITO_READS.out.bam  // channel: [ val(meta), ubam ]
        versions = ch_versions                // channel: [ versions.yml ]
}
