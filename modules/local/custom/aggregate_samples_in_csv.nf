process AGGREGATE_SAMPLES_IN_CSV {
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
