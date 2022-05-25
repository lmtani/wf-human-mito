include { BWA_ALIGN_FROM_UBAM        } from '../../modules/local/custom/bwa_align_from_ubam.nf'
include { FASTQ_TO_UBAM              } from '../../modules/local/picard/fastq_to_sam.nf'
include { PRINT_READS                } from '../../modules/local/gatk/print_reads.nf'
include { SELECT_MITO_READS          } from '../../modules/local/picard/revert_sam.nf'
include { SORT_SAM                   } from '../../modules/local/picard/sort_sam.nf'


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