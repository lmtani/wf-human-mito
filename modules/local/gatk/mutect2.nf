
process CALL_MUTECT {
    label "human_mito_mutect"
    clusterOptions "-C avx2"

    input:
        tuple val(sample_id), path(bam), path(bai)
        path mito_fasta
        path mito_dict
        path mito_index
        val prefix
        val mutect_extra_args
    output:
        tuple val(sample_id), \
            path("${prefix}.${sample_id}.vcf.gz"), \
            path("${prefix}.${sample_id}.vcf.gz.tbi"), \
            path("${prefix}.${sample_id}.vcf.gz.stats")

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
