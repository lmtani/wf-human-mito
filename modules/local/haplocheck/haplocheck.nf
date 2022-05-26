process HAPLOCHECK {
    tag "$sample_id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::haplocheck=1.3.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplocheck:1.3.3--h4a94de4_0':
        'quay.io/biocontainers/haplocheck:1.3.3--h4a94de4_0' }"

    input:

    tuple val(sample_id), path(vcf), path(vcf_index)

    output:
    tuple val(sample_id), path("*.txt"), emit: txt
    tuple val(sample_id), path("*.html"), emit: html

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    haplocheck --raw --out $sample_id $vcf
    """
}
