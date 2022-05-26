process SELECT_MITO_READS {

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
        tuple val(sample_id), path(bam)
        path fasta
        path dict
        path index

    output:
        tuple val(sample_id), path("${sample_id}.mito.unaligned.bam")

    script:
    """
    picard RevertSam \
        INPUT=${bam} \
        OUTPUT_BY_READGROUP=false \
        OUTPUT=${sample_id}.mito.unaligned.bam \
        VALIDATION_STRINGENCY=LENIENT \
        ATTRIBUTE_TO_CLEAR=FT \
        ATTRIBUTE_TO_CLEAR=CO \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=false
    """
}
