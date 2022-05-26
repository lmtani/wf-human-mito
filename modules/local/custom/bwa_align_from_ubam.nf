process BWA_ALIGN_FROM_UBAM {
    container "taniguti/wf-human-mito:${params.wfversion}"

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
        bwa mem -K 100000000 -p -v 3 -t !{task.cpus} -Y !{fasta} /dev/stdin - 2> >(tee !{sample_id}.bwa.stderr.log >&2) | \
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
