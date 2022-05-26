process CREATE_JSON {

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

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
