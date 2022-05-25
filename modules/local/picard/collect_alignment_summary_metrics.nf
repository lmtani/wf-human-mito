
process COLLECT_ALIGNMENT_METRICS {
    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

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
