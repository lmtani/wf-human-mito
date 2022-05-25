
process CREATE_JSON {
    label "human_mito"
    input:
        tuple \
            val(sample_id),
            path(contam_metrics),
            path(algn_metrics),
            path(wgs_metrics),
            path(theoretical_sensitivity),
            path(dup_metrics)

    output:
        path "${sample_id}.summary.json"

    script:
    """
    prepare_json.py $dup_metrics $wgs_metrics $algn_metrics $contam_metrics
    mv summary.json ${sample_id}.summary.json
    """
}
