process MARK_DUPLICATES {

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

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
