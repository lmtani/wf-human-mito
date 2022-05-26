process AGGREGATE_SAMPLES_IN_CSV {

    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
        path json

    output:
        path "all_samples.csv"

    script:
    """
    prepare_all_samples_table.py  "$json"
    """
}
