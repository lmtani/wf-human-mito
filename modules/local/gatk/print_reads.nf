process PRINT_READS {

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
        tuple val(sample_id), path(whole_bam), path(whole_bai)
        path fasta
        path index
        path dict

    output:
        tuple val(sample_id), path("mito.bam")

    script:
    """
    BAM=`find -L ./ -name "*.bam" -or -name "*.cram"`  # point to BAM, not BAI, in the gatk

    gatk PrintReads \
        -R $fasta \
        -L "chrM" \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -I \$BAM \
        -O mito.bam
    """
}
