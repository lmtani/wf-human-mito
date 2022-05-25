process FASTQ_TO_UBAM {
    label "human_mito"

    input:
        tuple \
            val(sampleId), \
            path(reads)

    output:
        tuple val(sampleId), path("${sampleId}.unmaped.bam")

    shell:
    """
    fastq_to_ubam.sh ${reads[0]} ${reads[1]} ${sampleId}
    """
}
