#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion.
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

nextflow.enable.dsl = 2


include { 
    BWA_ALIGN as ALIGN_STANDARD_MITO;
    BWA_ALIGN as ALIGN_SHIFTED_MITO;
    COLLECT_WGS_METRICS;
    CREATE_TABLE;
    BWA_ALIGN_FROM_UBAM;
    SELECT_MITO_READS;
    GET_CONTAMINATION;
    FILTER as INITIAL_FILTER;
    LIFTOVER_AND_COMBINE_VCFS;
    CALL_MUTECT as CALL_MUTECT_STANDARD;
    CALL_MUTECT as CALL_MUTECT_SHIFTED;
    SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES;
    FASTQ_TO_UBAM;
    MERGE_STATS} from './lib/align_and_call'

def helpMessage(){
    log.info """
Human Mitochondrial Analysis Workflow'

Usage:
    nextflow run lmtani/wf-human-mito [options]

Script Options:
    --fastq        DIR     Path to FASTQ directory. Quote is required. Ex: "/path/to/fastqs/*_R{1,2}*.fastq.gz" (required)
    --reference    Dir     Path to reference (GRCh38). BWA index need to be in same directory (required)
    --out_dir      DIR     Path for output (default: $params.out_dir)
"""
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "human_mito"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


// entrypoint workflow
workflow {

    if (params.help) {
        helpMessage()
        exit 1
    }

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }

    // Human Reference
    human_fasta = file("${params.reference}", type:'file', checkIfExists:true)
    human_dict = file("${human_fasta.getParent()}/${human_fasta.baseName}.dict", type:'file', checkIfExists:true)
    human_index = file("${params.reference}.fai", type:'file', checkIfExists:true)
    human_amb = file("${params.reference}.64.amb", type:'file', checkIfExists:true)
    human_ann = file("${params.reference}.64.ann", type:'file', checkIfExists:true)
    human_bwt = file("${params.reference}.64.bwt", type:'file', checkIfExists:true)
    human_sa = file("${params.reference}.64.sa", type:'file', checkIfExists:true)
    human_pac = file("${params.reference}.64.pac", type:'file', checkIfExists:true)
    human_alt = file("${params.reference}.64.alt", type:'file', checkIfExists:true)


    // Reference files - stored in github
    blacklist = file("$baseDir/data/blacklist_sites.hg38.chrM.bed", type:'file', checkIfExists:true)
    blacklist_index = file("$baseDir/data/blacklist_sites.hg38.chrM.bed.idx", type:'file', checkIfExists:true)
    shift_back_chain = file("$baseDir/data/ShiftBack.chain", type:'file', checkIfExists:true)

    mito_fasta = file("$baseDir/data/Homo_sapiens_assembly38.chrM.fasta", type:'file', checkIfExists:true)
    mito_dict = file("$baseDir/data/Homo_sapiens_assembly38.chrM.dict", type:'file', checkIfExists:true)
    mito_index = file("$baseDir/data/Homo_sapiens_assembly38.chrM.fasta.fai", type:'file', checkIfExists:true)
    mito_amb = file("$baseDir/data/Homo_sapiens_assembly38.chrM.fasta.amb", type:'file', checkIfExists:true)
    mito_ann = file("$baseDir/data/Homo_sapiens_assembly38.chrM.fasta.ann", type:'file', checkIfExists:true)
    mito_bwt = file("$baseDir/data/Homo_sapiens_assembly38.chrM.fasta.bwt", type:'file', checkIfExists:true)
    mito_sa = file("$baseDir/data/Homo_sapiens_assembly38.chrM.fasta.sa", type:'file', checkIfExists:true)
    mito_pac = file("$baseDir/data/Homo_sapiens_assembly38.chrM.fasta.pac", type:'file', checkIfExists:true)

    shifted_fasta = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta", type:'file', checkIfExists:true)
    shifted_dict = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.dict", type:'file', checkIfExists:true)
    shifted_index = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.fai", type:'file', checkIfExists:true)
    shifted_amb = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.amb", type:'file', checkIfExists:true)
    shifted_ann = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.ann", type:'file', checkIfExists:true)
    shifted_bwt = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.bwt", type:'file', checkIfExists:true)
    shifted_sa = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.sa", type:'file', checkIfExists:true)
    shifted_pac = file("$baseDir/data/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta.pac", type:'file', checkIfExists:true)

    reads = Channel.fromFilePairs("${params.fastq}", glob: true)

    FASTQ_TO_UBAM(reads)

    BWA_ALIGN_FROM_UBAM(
        FASTQ_TO_UBAM.out.sample_id,
        FASTQ_TO_UBAM.out.ubam,
        human_fasta,
        human_dict,
        human_index,
        human_amb,
        human_ann,
        human_bwt,
        human_pac,
        human_sa,
        human_alt
    )

    SELECT_MITO_READS(
        FASTQ_TO_UBAM.out.sample_id,
        BWA_ALIGN_FROM_UBAM.out.bam,
        BWA_ALIGN_FROM_UBAM.out.bai,
        human_fasta,
        human_dict,
        human_index,
    )

    ALIGN_STANDARD_MITO(
        SELECT_MITO_READS.out.ubam,
        FASTQ_TO_UBAM.out.sample_id,
        mito_fasta,
        mito_dict,
        mito_index,
        mito_amb,
        mito_ann,
        mito_bwt,
        mito_pac,
        mito_sa
    )

    ALIGN_SHIFTED_MITO(
        SELECT_MITO_READS.out.ubam,
        FASTQ_TO_UBAM.out.sample_id,
        shifted_fasta,
        shifted_dict,
        shifted_index,
        shifted_amb,
        shifted_ann,
        shifted_bwt,
        shifted_pac,
        shifted_sa
    )

    COLLECT_WGS_METRICS(
        ALIGN_STANDARD_MITO.out.bam,
        ALIGN_STANDARD_MITO.out.bai,
        mito_fasta,
        FASTQ_TO_UBAM.out.sample_id,
        300
    )


    CALL_MUTECT_STANDARD(
        ALIGN_STANDARD_MITO.out.bam,
        ALIGN_STANDARD_MITO.out.bai,
        mito_fasta,
        mito_dict,
        mito_index,
        "standard",
        FASTQ_TO_UBAM.out.sample_id,
        " -L chrM:576-16024 "
    )

    CALL_MUTECT_SHIFTED(
        ALIGN_SHIFTED_MITO.out.bam,
        ALIGN_SHIFTED_MITO.out.bai,
        shifted_fasta,
        shifted_dict,
        shifted_index,
        "shifted",
        FASTQ_TO_UBAM.out.sample_id,
        " -L chrM:8025-9144 "  
    )

    LIFTOVER_AND_COMBINE_VCFS(
        CALL_MUTECT_SHIFTED.out.vcf,
        CALL_MUTECT_STANDARD.out.vcf,
        FASTQ_TO_UBAM.out.sample_id,
        mito_fasta,
        mito_index,
        mito_dict,
        shift_back_chain
    )

    MERGE_STATS(
        CALL_MUTECT_SHIFTED.out.stats,
        CALL_MUTECT_STANDARD.out.stats
    )

    INITIAL_FILTER(
        mito_fasta,
        mito_index,
        mito_dict,
        LIFTOVER_AND_COMBINE_VCFS.out.merged_vcf,
        LIFTOVER_AND_COMBINE_VCFS.out.merged_vcf_idx,
        MERGE_STATS.out,
        FASTQ_TO_UBAM.out.sample_id,
        blacklist,
        blacklist_index
    )

    SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES(
        mito_fasta,
        mito_index,
        mito_dict,
        INITIAL_FILTER.out.vcf,
        INITIAL_FILTER.out.tbi,
        FASTQ_TO_UBAM.out.sample_id,
    )

    GET_CONTAMINATION(
        SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out.vcf
    )

    CREATE_TABLE(
        ALIGN_STANDARD_MITO.out.metrics,
        COLLECT_WGS_METRICS.out.metrics,
        GET_CONTAMINATION.out.contamination_file,
        FASTQ_TO_UBAM.out.sample_id
    )

    output(
        ALIGN_STANDARD_MITO.out.bam.concat(
            ALIGN_STANDARD_MITO.out.bai, 
            // ALIGN_STANDARD_MITO.out.metrics,
            // COLLECT_WGS_METRICS.out.sensitivity,
            // COLLECT_WGS_METRICS.out.metrics,
            LIFTOVER_AND_COMBINE_VCFS.out.rejected,
            // SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out.vcf,
            // SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out.tbi,
            INITIAL_FILTER.out.vcf,
            INITIAL_FILTER.out.tbi,
            // MERGE_STATS.out,
            // GET_CONTAMINATION.out.major_hg,
            // GET_CONTAMINATION.out.minor_hg,
            // LIFTOVER_AND_COMBINE_VCFS.out.merged_vcf,
            CREATE_TABLE.out,
            // GET_CONTAMINATION.out.major_level,
            // GET_CONTAMINATION.out.minor_level,
            // GET_CONTAMINATION.out.contamination_file
        )
    )
}
