process MERGE_VCFS {

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
        tuple val(meta), path(vcf1), path(vcf2)

    output:
        tuple val(meta), path("${meta.id}.merged.vcf")       , emit: vcf
        tuple val(meta), path("${meta.id}.merged.vcf.idx")   , emit: idx
        path "versions.yml"                                  , emit: versions

    script:
    """
    picard MergeVcfs \
        I=${vcf1} \
        I=${vcf2} \
        O=${meta.id}.merged.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MergeVcfs --version 2> >(grep -v LC_ALL))
    END_VERSIONS
    """
}
