process PRINT_READS {
    label "human_mito"

    input:
        tuple val(sample_id), path(whole_bam), path(whole_bai)
        path fasta
        path index
        path dict

    output:
        tuple val(sample_id), path("mito.bam")

    script:
    """
    gatk PrintReads \
        -R $fasta \
        -L "chrM" \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -I $whole_bam \
        -O mito.bam
    """
}
