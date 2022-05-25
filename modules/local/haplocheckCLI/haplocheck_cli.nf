process GET_CONTAMINATION {
    container "us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124"  //TODO: change to pure haplocheck

    input:
        tuple val(sample_id), path(vcf), path(vcf_index)
    output:
        tuple val(sample_id), path("output-noquotes")

    script:
    """
    get_contamination.sh $vcf
    """
}
