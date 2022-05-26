process LIFTOVER_VCF {
    // label "human_mito"

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
        tuple val(sample_id), path(vcf), path(tbi)
        path mito_fasta
        path mito_fasta_index
        path mito_dict
        path chain

    output:
        tuple val(sample_id), path("${sample_id}.lifted.vcf.gz"), path("${sample_id}.lifted.vcf.gz.tbi")

    script:
    """
    picard LiftoverVcf \
        I=${vcf} \
        O=${sample_id}.lifted.vcf.gz \
        R=${mito_fasta} \
        CHAIN=${chain} \
        REJECT=${sample_id}.rejected.vcf
    """
}
