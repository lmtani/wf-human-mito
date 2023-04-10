process SELECT_MITO_READS {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
        tuple val(meta), path(bam)
        path fasta
        path dict
        path index
        val restore_hardclips

    output:
        tuple val(meta), path("${meta.id}.mito.unaligned.bam"), emit: bam
        path "versions.yml"                                   , emit: versions

    script:
    if (!restore_hardclips) {
        args = "RESTORE_HARDCLIPS=false"
    } else {
        args = ""
    }
    """
    picard RevertSam \
        INPUT=${bam} \
        OUTPUT_BY_READGROUP=false \
        OUTPUT=${meta.id}.mito.unaligned.bam \
        VALIDATION_STRINGENCY=LENIENT \
        ATTRIBUTE_TO_CLEAR=FT \
        ATTRIBUTE_TO_CLEAR=CO \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=false \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard RevertSam --version 2> >(grep -v LC_ALL))
    END_VERSIONS
    """
}
