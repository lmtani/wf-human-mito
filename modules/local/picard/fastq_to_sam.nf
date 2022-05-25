process FASTQ_TO_UBAM {
    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

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
