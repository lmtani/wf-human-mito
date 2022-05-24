process SORT_SAM {
    label "human_mito"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bai")

    script:
    """
    picard SortSam \
        INPUT="${bam}" \
        OUTPUT="${sample_id}.sorted.bam" \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        MAX_RECORDS_IN_RAM=300000
    """
}
