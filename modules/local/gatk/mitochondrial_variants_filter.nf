process FILTER_MUTECT_CALLS {
    label "human_mito"
    input:
        tuple val(sample_id), \
            path(merged_vcf), \
            path(merged_vcf_index), \
            path(combined_stats)
        path mito_fasta
        path mito_index
        path mito_dict
        path blacklisted_sites
        path blacklisted_sites_index

    output:
        tuple val(sample_id), \
            path("${sample_id}.vcf.gz"), \
            path("${sample_id}.vcf.gz.tbi")

    """
    gatk FilterMutectCalls \
        -V ${merged_vcf} \
        -R ${mito_fasta} \
        -O filtered.vcf \
        --stats ${combined_stats} \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        --min-allele-fraction 0

    gatk VariantFiltration -V filtered.vcf \
        -O ${sample_id}.vcf.gz \
        --apply-allele-specific-filters \
        --mask ${blacklisted_sites} \
        --mask-name "blacklisted_site"
    """
}
