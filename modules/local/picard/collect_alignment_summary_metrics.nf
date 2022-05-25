
process COLLECT_ALIGNMENT_METRICS {
    label "human_mito"

    input:
        tuple val(sample_id), path(bam), path(bai)
        path reference
        path reference_dict
        path reference_index

    output:
        tuple val(sample_id), path("${sample_id}.algn_metrics.txt")

    script:
    """
    picard CollectAlignmentSummaryMetrics \
        R=$reference \
        I=$bam \
        O=${sample_id}.algn_metrics.txt
    """
}
