process MERGE_VCFS {
    label "human_mito"
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
