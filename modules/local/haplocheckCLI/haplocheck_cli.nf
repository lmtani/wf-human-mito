process GET_CONTAMINATION {
    label "mtdnaserver"
    input:
        tuple val(sample_id), path(vcf), path(vcf_index)
    output:
        tuple val(sample_id), path("output-noquotes")

    script:
    """
    get_contamination.sh $vcf
    """
}
