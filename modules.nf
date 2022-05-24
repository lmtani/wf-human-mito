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
        tuple val(sample_id), path("${sample_id}.temp.bam")
    
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
    """
}

// process BWA_ALIGN {
//     label "human_mito"

//     input:
//         tuple val(sampleId), path(input_bam)
//         path mito_fasta
//         path mito_dict
//         path mito_index
//         path mito_amb
//         path mito_ann
//         path mito_bwt
//         path mito_pac
//         path mito_sa
//     output:
//         tuple \
//             val(sampleId), \
//             path("${sampleId}.bam"), \
//             path("${sampleId}.bai"), \
//             path("${sampleId}.dup.metrics")

//     shell:
//     """
//     picard SamToFastq \
//         INPUT=!{input_bam} \
//         FASTQ=/dev/stdout \
//         INTERLEAVE=true \
//         NON_PF=true | \
//     bwa mem -K 100000000 -p -v 3 -t 2 -Y !{mito_fasta} /dev/stdin - 2> >(tee !{sampleId}.bwa.stderr.log >&2) | \
//     picard MergeBamAlignment \
//         VALIDATION_STRINGENCY=SILENT \
//         EXPECTED_ORIENTATIONS=FR \
//         ATTRIBUTES_TO_RETAIN=X0 \
//         ATTRIBUTES_TO_REMOVE=NM \
//         ATTRIBUTES_TO_REMOVE=MD \
//         ALIGNED_BAM=/dev/stdin \
//         UNMAPPED_BAM=!{input_bam} \
//         OUTPUT=mba.bam \
//         REFERENCE_SEQUENCE=!{mito_fasta} \
//         PAIRED_RUN=true \
//         SORT_ORDER="unsorted" \
//         IS_BISULFITE_SEQUENCE=false \
//         ALIGNED_READS_ONLY=false \
//         CLIP_ADAPTERS=false \
//         MAX_RECORDS_IN_RAM=2000000 \
//         ADD_MATE_CIGAR=true \
//         MAX_INSERTIONS_OR_DELETIONS=-1 \
//         PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
//         UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
//         ALIGNER_PROPER_PAIR_FLAGS=true \
//         UNMAP_CONTAMINANT_READS=true \
//         ADD_PG_TAG_TO_READS=false

//     picard MarkDuplicates \
//         INPUT=mba.bam \
//         OUTPUT=md.bam \
//         METRICS_FILE=!{sampleId}.dup.metrics \
//         VALIDATION_STRINGENCY=SILENT \
//         OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
//         ASSUME_SORT_ORDER="queryname" \
//         CLEAR_DT="false" \
//         ADD_PG_TAG_TO_READS=false

//     picard SortSam \
//         INPUT=md.bam \
//         OUTPUT=!{sampleId}.bam \
//         SORT_ORDER="coordinate" \
//         CREATE_INDEX=true \
//         MAX_RECORDS_IN_RAM=300000
//     """
// }

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
        tuple val(sample_id), path("combined.stats")
    
    script:
    """
    gatk MergeMutectStats \
        --stats ${shifted_stats} \
        --stats ${stantard_stats} \
        -O combined.stats
    """
}

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
        tuple val(sample_id), \
            path("${sample_id}.rejected.vcf"), \
            path("${sample_id}.merged.vcf"), \
            path("${sample_id}.merged.vcf.idx")

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
        tuple val(sample_id), \
            path(rejected_vcf), \
            path(merged_vcf), \
            path(merged_vcf_index), \
            path(combined_stats)
        path mito_fasta
        path mito_index
        path mito_dict
        path blacklisted_sites
        path blacklisted_sites_index

    output:
        tuple val(sample_id), \
            path("${sample_id}.vcf.gz"), \
            path("${sample_id}.vcf.gz.tbi")

    """
    gatk FilterMutectCalls \
        -V ${merged_vcf} \
        -R ${mito_fasta} \
        -O filtered.vcf \
        --stats ${combined_stats} \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        --min-allele-fraction 0

    gatk VariantFiltration -V filtered.vcf \
        -O ${sample_id}.vcf.gz \
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
        tuple val(sample_id), path(filtered_vcf), path(filtered_vcf_index)

    output:
        tuple val(sample_id), \
            path("${sample_id}.pass.vcf.gz"), \
            path("${sample_id}.pass.vcf.gz.tbi")

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
        tuple val(sample_id), path(vcf), path(vcf_index)
    output:
        tuple val(sample_id), path("output-noquotes")

    script:
    """
    get_contamination.sh $vcf
    """
}

process CREATE_JSON {
    label "human_mito"
    input:
        tuple val(sample_id), \
            path(vcf), \
            path(vcf_idx), \
            path(contam_metrics), \
            path(bam), \
            path(bai), \
            path(dup_metrics), \
            path(algn_metrics), \
            path(theoretical_sensitivity), \
            path(wgs_metrics)

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
