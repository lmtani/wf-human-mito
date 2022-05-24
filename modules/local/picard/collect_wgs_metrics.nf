process COLLECT_WGS_METRICS {
    label "human_mito"

    input:
        tuple val(sample_id), path(bam), path(bai)
        path reference
        val readLen
    output:
        tuple val(sample_id), \
            path("${sample_id}.theoretical_sensitivity.txt"), \
            path("${sample_id}.metrics.txt")

    script:
    """
    picard CollectWgsMetrics \
        INPUT=${bam} \
        VALIDATION_STRINGENCY=SILENT \
        REFERENCE_SEQUENCE=${reference} \
        OUTPUT=${sample_id}.metrics.txt \
        USE_FAST_ALGORITHM=true \
        READ_LENGTH=${readLen} \
        INCLUDE_BQ_HISTOGRAM=true \
        COVERAGE_CAP=100000 \
        THEORETICAL_SENSITIVITY_OUTPUT=${sample_id}.theoretical_sensitivity.txt
    """
}