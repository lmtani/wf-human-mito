process SORT_SAM {
    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

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
