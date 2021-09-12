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


process SELECT_MITO_READS {
    label "human_mito"

    input:
        tuple val(sample_id), path(whole_bam), path(whole_bai)
        path fasta
        path dict
        path index

    output:
        tuple val(sample_id), path("${sample_id}.mito.unaligned.bam")

    script:
    """
    gatk PrintReads \
        -R $fasta \
        -L "chrM" \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -I $whole_bam \
        -O mito.bam

    picard RevertSam \
        INPUT=mito.bam \
        OUTPUT_BY_READGROUP=false \
        OUTPUT=${sample_id}.mito.unaligned.bam \
        VALIDATION_STRINGENCY=LENIENT \
        ATTRIBUTE_TO_CLEAR=FT \
        ATTRIBUTE_TO_CLEAR=CO \
        SORT_ORDER=queryname \
        RESTORE_ORIGINAL_QUALITIES=false
    """
}


process BWA_ALIGN_FROM_UBAM {
    label "human_mito"

    input:
        tuple val(sample_id), path(ubam)
        path fasta
        path dict
        path index
        path amb
        path ann
        path bwt
        path pac
        path sa
        path ref_alt
    output:
        tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bai")
    
    shell:
    """
    picard SamToFastq \
            INPUT=!{ubam} \
            FASTQ=/dev/stdout \
            INTERLEAVE=true \
            NON_PF=true | \
        bwa mem -K 100000000 -p -v 3 -t 5 -Y !{fasta} /dev/stdin - 2> >(tee !{sample_id}.bwa.stderr.log >&2) | \
        picard MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM=!{ubam} \
            OUTPUT=!{sample_id}.temp.bam \
            REFERENCE_SEQUENCE=!{fasta} \
            SORT_ORDER="unsorted" \
            IS_BISULFITE_SEQUENCE=false \
            ALIGNED_READS_ONLY=false \
            CLIP_ADAPTERS=false \
            MAX_RECORDS_IN_RAM=2000000 \
            ADD_MATE_CIGAR=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=true \
            ADD_PG_TAG_TO_READS=false

    picard SortSam \
        INPUT="!{sample_id}.temp.bam" \
        OUTPUT="!{sample_id}.sorted.bam" \
        SORT_ORDER="coordinate" \
        CREATE_INDEX=true \
        MAX_RECORDS_IN_RAM=300000
    """
}

process BWA_ALIGN {
    label "human_mito"

    input:
        tuple val(sampleId), path(input_bam)
        path mito_fasta
        path mito_dict
        path mito_index
        path mito_amb
        path mito_ann
        path mito_bwt
        path mito_pac
        path mito_sa
    output:
        tuple \
            val(sampleId), \
            path("${sampleId}.bam"), \
            path("${sampleId}.bai"), \
            path("${sampleId}.dup.metrics")


    
    shell:
    """
    picard SamToFastq \
      INPUT=!{input_bam} \
      FASTQ=/dev/stdout \
      INTERLEAVE=true \
      NON_PF=true | \
    bwa mem -K 100000000 -p -v 3 -t 2 -Y !{mito_fasta} /dev/stdin - 2> >(tee !{sampleId}.bwa.stderr.log >&2) | \
    picard MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ATTRIBUTES_TO_REMOVE=NM \
      ATTRIBUTES_TO_REMOVE=MD \
      ALIGNED_BAM=/dev/stdin \
      UNMAPPED_BAM=!{input_bam} \
      OUTPUT=mba.bam \
      REFERENCE_SEQUENCE=!{mito_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=true \
      ADD_PG_TAG_TO_READS=false

    picard MarkDuplicates \
      INPUT=mba.bam \
      OUTPUT=md.bam \
      METRICS_FILE=!{sampleId}.dup.metrics \
      VALIDATION_STRINGENCY=SILENT \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false

    picard SortSam \
      INPUT=md.bam \
      OUTPUT=!{sampleId}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      MAX_RECORDS_IN_RAM=300000
    """
}

process COLLECT_ALIGNMENT_METRICS {
    label "human_mito"

    input:
        tuple val(sample_id), path(bam), path(bai), path(dup_stats)
        path reference
        path reference_dict
        path reference_index
    
    output:
        path "algn_metrics.txt"

    script:
    """
    picard CollectAlignmentSummaryMetrics \
          R=$reference \
          I=$bam \
          O=algn_metrics.txt
    """
}

process COLLECT_WGS_METRICS {
    label "human_mito"

    input:
        tuple val(sample_id), path(bam), path(bai), path(dup_stats)
        path reference
        val readLen
    output:
        path "${sample_id}.theoretical_sensitivity.txt", emit: sensitivity
        path "${sample_id}.metrics.txt", emit: metrics
    script:
    """
    picard CollectWgsMetrics \
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
    clusterOptions "-C avx2"

    input:
        tuple val(sample_id), path(bam), path(bai), path(dup_stats)
        path mito_fasta
        path mito_dict
        path mito_index
        val prefix
        val mutect_extra_args
    output:
        tuple val(sample_id), \
            path("${prefix}.${sample_id}.vcf.gz"), \
            path("${prefix}.${sample_id}.vcf.gz.tbi"), \
            path("${prefix}.${sample_id}.vcf.gz.stats")

    script:
    """
    gatk Mutect2 \
        -R ${mito_fasta} \
        -I ${bam} \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O "${prefix}.${sample_id}.vcf.gz" \
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
        tuple \
            val(sample_id), \
            path(standard_vcf), \
            path(standard_tbi), \
            path(stantard_stats), \
            path(shifted_vcf), \
            path(shifted_tbi), \
            path(shifted_stats)

    output:
        path "combined.stats"
    
    script:
    """
    gatk MergeMutectStats \
        --stats ${shifted_stats} \
        --stats ${stantard_stats} \
        -O combined.stats
    """
}
// [GEF495-001, 
// /home/taniguti/github/wf-human-mito/work/40/441c17ead6aa98e50d7b80cda0b003/standard.GEF495-001.vcf.gz, 
// /home/taniguti/github/wf-human-mito/work/40/441c17ead6aa98e50d7b80cda0b003/standard.GEF495-001.vcf.gz.tbi, 
// /home/taniguti/github/wf-human-mito/work/40/441c17ead6aa98e50d7b80cda0b003/standard.GEF495-001.vcf.gz.stats, 
// /home/taniguti/github/wf-human-mito/work/05/f7ef2ad5f1edc682bbd8d5c49aa290/shifted.GEF495-001.vcf.gz, 
// /home/taniguti/github/wf-human-mito/work/05/f7ef2ad5f1edc682bbd8d5c49aa290/shifted.GEF495-001.vcf.gz.tbi, 
// /home/taniguti/github/wf-human-mito/work/05/f7ef2ad5f1edc682bbd8d5c49aa290/shifted.GEF495-001.vcf.gz.stats]

process LIFTOVER_AND_COMBINE_VCFS {
    label "human_mito"
    input:
        tuple \
            val(sample_id), \
            path(standard_vcf), \
            path(standard_tbi), \
            path(stantard_stats), \
            path(shifted_vcf), \
            path(shifted_tbi), \
            path(shifted_stats)
        path mito_fasta
        path mito_fasta_index
        path mito_dict
        path shift_back_chain

    output:
        path "${sample_id}.rejected.vcf", emit: rejected
        path "${sample_id}.merged.vcf", emit: merged_vcf
        path "${sample_id}.merged.vcf.idx", emit: merged_vcf_idx
        val sample_id, emit: basename
    script:
    """
    picard LiftoverVcf \
      I=${shifted_vcf} \
      O=${sample_id}.shifted_back.vcf \
      R=${mito_fasta} \
      CHAIN=${shift_back_chain} \
      REJECT=${sample_id}.rejected.vcf

    picard MergeVcfs \
      I=${sample_id}.shifted_back.vcf \
      I=${standard_vcf} \
      O=${sample_id}.merged.vcf
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
        val basename, emit: basename

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

process CREATE_JSON {
    label "human_mito"
    input:
        tuple val(sample_id), path(bam), path(bai), path(dup_metrics)
        path wgs_metrics
        path contam_metrics
        path algn_metrics

    output:
        path "${sample_id}.summary.json"

    script:
    """
    prepare_json.py $dup_metrics $wgs_metrics $algn_metrics $contam_metrics
    mv summary.json ${sample_id}.summary.json
    """
}


process CREATE_ALL_SAMPLES_CSV {
    label "human_mito"
    input:
        path json

    output:
        path "all_samples.csv"

    script:
    """
    prepare_all_samples_table.py  "$json"
    """
}
