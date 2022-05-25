
process SELECT_VARIANTS {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
        tuple val(sample_id), path(vcf), path(tbi)

    output:
        tuple val(sample_id), \
            path("${sample_id}.filtered.vcf.gz"), \
            path("${sample_id}.filtered.vcf.gz.tbi")

    script:
    """
    gatk SelectVariants \
        -V ${vcf} \
        -O ${sample_id}.filtered.vcf.gz \
        --exclude-filtered
    """
}
