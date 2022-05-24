process MARK_DUPLICATES {
    label "human_mito"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("${sample_id}.markdup.bam"), emit: bam
        tuple val(sample_id), path("${sample_id}.dup.metrics"), emit: metrics

    script:
    """
    picard MarkDuplicates \
        INPUT=${bam} \
        OUTPUT=${sample_id}.markdup.bam \
        METRICS_FILE=${sample_id}.dup.metrics \
        VALIDATION_STRINGENCY=SILENT \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER="queryname" \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false
    """
}
