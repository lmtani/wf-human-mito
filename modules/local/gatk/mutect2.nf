
process CALL_MUTECT {
    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    clusterOptions "-C avx2"  // TODO: improve visibility of this requirement

    input:
        tuple val(sample_id), path(bam), path(bai)
        path mito_fasta
        path mito_dict
        path mito_index
        val prefix
        val mutect_extra_args
    output:
        tuple val(sample_id), path("${prefix}.${sample_id}.vcf.gz"), path("${prefix}.${sample_id}.vcf.gz.tbi"), emit: vcf
        tuple val(sample_id), path("${prefix}.${sample_id}.vcf.gz.stats"), emit: stats

    script:
    """
    gatk Mutect2 \
        -R ${mito_fasta} \
        -I ${bam} \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O "${prefix}.${sample_id}.vcf.gz" \
        ${mutect_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0
    """
}
