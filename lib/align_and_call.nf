process BWA_ALIGN {
    label "human_mito"

    input:
        tuple \
            val(sampleId), \
            path(reads)
        path mito_fasta
        path mito_dict
        path mito_index
        path mito_amb
        path mito_ann
        path mito_bwt
        path mito_pac
        path mito_sa
    output:
        val "${sampleId}", emit: sample_id
        path "${sampleId}.sorted.bam", emit: bam
        path "${sampleId}.sorted.bai", emit: bai
        path "${sampleId}.dup.metrics", emit: metrics

    
    script:
    """
    alignment_pipeline.sh ${mito_fasta} ${reads[0]} ${reads[1]} ${sampleId}
    """
}

process COLLECT_WGS_METRICS {
    label "human_mito"

    input:
        path bam
        path bai
        path reference
        val sample_id
        val readLen
    output:
        path "${sample_id}.theoretical_sensitivity.txt", emit: sensitivity
        path "${sample_id}.metrics.txt", emit: metrics
    script:
    """
    picard \
        CollectWgsMetrics \
        INPUT=${bam} \
        VALIDATION_STRINGENCY=SILENT \
        REFERENCE_SEQUENCE=${reference} \
        OUTPUT=${sample_id}.metrics.txt \
        USE_FAST_ALGORITHM=true \
        READ_LENGTH=${readLen} \
        INCLUDE_BQ_HISTOGRAM=true \
        COVERAGE_CAP=100000 \
        THEORETICAL_SENSITIVITY_OUTPUT=${sample_id}.theoretical_sensitivity.txt
    """
}

process CALL_MUTECT {
    label "human_mito"

    input:
        path bam
        path bai
        path mito_fasta
        path mito_dict
        path mito_index
        val prefix
        val basename
        val mutect_extra_args
    output:
        path "${prefix}.${basename}.vcf.gz", emit: vcf
        path "${prefix}.${basename}.vcf.gz.tbi", emit: tbi
        path "${prefix}.${basename}.vcf.gz.stats", emit: stats
    script:
    """
    gatk Mutect2 \
        -R ${mito_fasta} \
        -I ${bam} \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O "${prefix}.${basename}.vcf.gz" \
        ${mutect_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0
    """
}


process MERGE_STATS {
    label "human_mito"

    input:
        path shifted_stats
        path standard_stats

    output:
        path "combined.stats"
    
    script:
    """
    gatk MergeMutectStats \
        --stats ${shifted_stats} \
        --stats ${standard_stats} \
        -O combined.stats
    """
}


process LIFTOVER_AND_COMBINE_VCFS {
    label "human_mito"
    input:
        path shifted_vcf
        path vcf
        val basename
        path mito_fasta
        path mito_fasta_index
        path mito_dict
        path shift_back_chain

    output:
        path "${basename}.rejected.vcf", emit: rejected
        path "${basename}.merged.vcf", emit: merged_vcf
        path "${basename}.merged.vcf.idx", emit: merged_vcf_idx
    script:
    """
    picard LiftoverVcf \
      I=${shifted_vcf} \
      O=${basename}.shifted_back.vcf \
      R=${mito_fasta} \
      CHAIN=${shift_back_chain} \
      REJECT=${basename}.rejected.vcf

    picard MergeVcfs \
      I=${basename}.shifted_back.vcf \
      I=${vcf} \
      O=${basename}.merged.vcf
    """
}

process FILTER {
    label "human_mito"
    input:
        path mito_fasta
        path mito_index
        path mito_dict
        path raw_vcf
        path raw_vcf_index
        path raw_vcf_stats
        val basename
        path blacklisted_sites
        path blacklisted_sites_index

    output:
        path "${basename}.vcf.gz", emit: vcf
        path "${basename}.vcf.gz.tbi", emit: tbi

    """
    gatk FilterMutectCalls \
        -V ${raw_vcf} \
        -R ${mito_fasta} \
        -O filtered.vcf \
        --stats ${raw_vcf_stats} \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        --min-allele-fraction 0

    gatk VariantFiltration -V filtered.vcf \
        -O ${basename}.vcf.gz \
        --apply-allele-specific-filters \
        --mask ${blacklisted_sites} \
        --mask-name "blacklisted_site"
    """
}

process SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES {
    label "human_mito"
    input:
        path mito_fasta
        path mito_index
        path mito_dict
        path filtered_vcf
        path filtered_vcf_index
        val sample_id

    output:
        path "${sample_id}.pass.vcf.gz", emit: vcf
        path "${sample_id}.pass.vcf.gz.tbi", emit: tbi

    script:
    """
    gatk LeftAlignAndTrimVariants \
      -R ${mito_fasta} \
      -V ${filtered_vcf} \
      -O split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

      gatk SelectVariants \
        -V split.vcf \
        -O ${sample_id}.pass.vcf.gz \
        --exclude-filtered
    """
}

process GET_CONTAMINATION {
    label "mtdnaserver"
    input:
        path vcf
    output:
        path "output-noquotes", emit: contamination_file
        path "contamination.txt", emit: has_contamination
        path "major_hg.txt", emit: major_hg
        path "minor_hg.txt", emit: minor_hg
        path "mean_het_major.txt", emit: major_level
        path "mean_het_minor.txt", emit: minor_level

    script:
    """
    get_contamination.sh $vcf
    """
}

process CREATE_TABLE {
    label "human_mito"
    input:
        file minor_hg
        file major_hg

    output:
        path "teste"

    script:
    minor = minor_hg.getText()
    """
    echo "$minor_hg,$major_hg" > teste
    """
}