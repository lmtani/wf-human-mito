process CREATE_JSON {

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
        tuple val(meta), path(contam_metrics), path(dup_metrics)

    output:
        path "${meta.id}.summary.json"

    script:
    """
    prepare_json.py $dup_metrics $contam_metrics
    mv summary.json ${meta.id}.summary.json
    """
}
