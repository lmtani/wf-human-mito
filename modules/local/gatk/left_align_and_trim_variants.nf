process LEFT_ALIGN_AND_TRIM_VARIANTS {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
        path mito_fasta
        path mito_index
        path mito_dict
        tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_index)

    output:
        tuple val(sample_id), \
            path("${sample_id}.split.vcf.gz"), \
            path("${sample_id}.split.vcf.gz.tbi")

    script:
    """
    gatk LeftAlignAndTrimVariants \
        -R ${mito_fasta} \
        -V ${filtered_vcf} \
        -O ${sample_id}.split.vcf.gz \
        --split-multi-allelics \
        --dont-trim-alleles \
        --keep-original-ac
    """
}
