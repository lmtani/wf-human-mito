process MERGE_STATS {
    label "human_mito"

    input:
        tuple val(sample_id), path(stantard_stats), path(shifted_stats)

    output:
        tuple val(sample_id), path("combined.stats")
    
    script:
    """
    gatk MergeMutectStats \
        --stats ${shifted_stats} \
        --stats ${stantard_stats} \
        -O combined.stats
    """
}