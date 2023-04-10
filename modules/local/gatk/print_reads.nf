process PRINT_READS {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
        tuple val(meta), path(whole_bam), path(whole_bai)
        path fasta
        path index
        path dict

    output:
        tuple val(meta), path("mito.bam") , emit: bam
        path "versions.yml"                  , emit: versions
    script:
    """
    BAM=`find -L ./ -name "*.bam" -or -name "*.cram"`  # point to BAM, not BAI, in the gatk

    # Uses the ref.dict to infer if the mitochondrial genome is named MT or chrMT
    CHROM_NAME=`grep -q "MT" $dict && echo "MT" || echo "chrM"`

    gatk PrintReads \
        -R $fasta \
        -L "\$CHROM_NAME" \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -I \$BAM \
        -O mito.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
