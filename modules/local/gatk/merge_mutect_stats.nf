process MERGE_STATS {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
        tuple val(sample_id), path(stantard_stats), path(shifted_stats)

    output:
        tuple val(sample_id), path("combined.stats")

    script:
    """
    gatk MergeMutectStats \
        --stats ${shifted_stats} \
        --stats ${stantard_stats} \
        -O combined.stats
    """
}
