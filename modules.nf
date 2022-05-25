


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

process CREATE_JSON {
    label "human_mito"
    input:
        tuple val(sample_id), \
            path(vcf), \
            path(vcf_idx), \
            path(contam_metrics), \
            path(bam), \
            path(bai), \
            path(dup_metrics), \
            path(algn_metrics), \
            path(theoretical_sensitivity), \
            path(wgs_metrics)

    output:
        path "${sample_id}.summary.json"

    script:
    """
    prepare_json.py $dup_metrics $wgs_metrics $algn_metrics $contam_metrics
    mv summary.json ${sample_id}.summary.json
    """
}


process CREATE_ALL_SAMPLES_CSV {
    label "human_mito"
    input:
        path json

    output:
        path "all_samples.csv"

    script:
    """
    prepare_all_samples_table.py  "$json"
    """
}
