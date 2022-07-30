//
// Aligns paired FASTQ reads into the human reference and selects
// reads from mitochondrial genome.
//
include { BWA_ALIGN_FROM_UBAM        } from '../../modules/local/custom/bwa_align_from_ubam.nf'
include { GATK4_FASTQTOSAM           } from '../../modules/local/gatk/fastq_to_sam.nf'
include { PRINT_READS                } from '../../modules/local/gatk/print_reads.nf'
include { SELECT_MITO_READS          } from '../../modules/local/picard/revert_sam.nf'
include { SORT_SAM                   } from '../../modules/local/picard/sort_sam.nf'


workflow separate_mitochondrion {
    take:
        reads             // channel: [ val(sample_id), [ path(fq_r1), path(fq_r2) ] ]
        alignments        // channel: [ val(sample_id), [ path(bam), path(bai) ] ]
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

        GATK4_FASTQTOSAM(reads)
        BWA_ALIGN_FROM_UBAM(GATK4_FASTQTOSAM.out, fasta, dict, index, amb, ann, bwt, pac, sa, alt)
        SORT_SAM(BWA_ALIGN_FROM_UBAM.out)

        mito_reads_ch = SORT_SAM.out.mix(alignments)
        PRINT_READS(mito_reads_ch, fasta, index, dict)

        SELECT_MITO_READS(PRINT_READS.out, fasta, dict, index, restore_hardclips)

    emit:
        SELECT_MITO_READS.out  // channel: [ val(sample_id), ubam ]
}
