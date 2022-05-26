process MERGE_VCFS {
    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
        tuple val(sample_id), path(vcf1), path(tbi1), path(vcf2), path(tbi2)

    output:
        tuple val(sample_id), \
            path("${sample_id}.merged.vcf"), \
            path("${sample_id}.merged.vcf.idx")

    script:
    """
    picard MergeVcfs \
        I=${vcf1} \
        I=${vcf2} \
        O=${sample_id}.merged.vcf
    """
}
