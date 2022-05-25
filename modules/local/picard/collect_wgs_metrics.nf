process COLLECT_WGS_METRICS {
    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

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
