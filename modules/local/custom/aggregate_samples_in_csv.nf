process AGGREGATE_SAMPLES_IN_CSV {

    container "taniguti/wf-human-mito:${params.wfversion}"

    input:
        path json

    output:
        path "all_samples.csv"

    script:
    """
    prepare_all_samples_table.py  "$json"
    """
}
