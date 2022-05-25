
process SELECT_VARIANTS {
    label "human_mito"
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
