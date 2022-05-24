process SELECT_MITO_READS {
    label "human_mito"

    input:
        tuple val(sample_id), path(bam)
        path fasta
        path dict
        path index

    output:
        tuple val(sample_id), path("${sample_id}.mito.unaligned.bam")

    script:
    """
    picard RevertSam \
        INPUT=${bam} \
        OUTPUT_BY_READGROUP=false \
        OUTPUT=${sample_id}.mito.unaligned.bam \
        VALIDATION_STRINGENCY=LENIENT \
        ATTRIBUTE_TO_CLEAR=FT \
        ATTRIBUTE_TO_CLEAR=CO \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=false
    """
}
