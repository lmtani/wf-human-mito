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
    GET_CONTAMINATION;
    FILTER as INITIAL_FILTER;
    LIFTOVER_AND_COMBINE_VCFS;
    CALL_MUTECT as CALL_MUTECT_STANDARD;
    CALL_MUTECT as CALL_MUTECT_SHIFTED;
    SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES;
    MERGE_STATS} from './lib/align_and_call'

def helpMessage(){
    log.info """
Human Mitochondrial Analysis Workflow'

Usage:
    nextflow run lmtani/wf-human-mito [options]

Script Options:
    --fastq        DIR     Path to FASTQ directory. Ex: "/path/to/fastqs/*_R{1,2}*.fastq.gz" (required)
    --out_dir      DIR     Path for output (default: $params.out_dir)
"""
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfcas9"
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

    ALIGN_STANDARD_MITO(
        reads,
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
        reads,
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
        ALIGN_STANDARD_MITO.out.sample_id,
        300
    )


    CALL_MUTECT_STANDARD(
        ALIGN_STANDARD_MITO.out.bam,
        ALIGN_STANDARD_MITO.out.bai,
        mito_fasta,
        mito_dict,
        mito_index,
        "standard",
        ALIGN_STANDARD_MITO.out.sample_id,
        " -L chrM:576-16024 "
    )

    CALL_MUTECT_SHIFTED(
        ALIGN_SHIFTED_MITO.out.bam,
        ALIGN_SHIFTED_MITO.out.bai,
        shifted_fasta,
        shifted_dict,
        shifted_index,
        "shifted",
        ALIGN_SHIFTED_MITO.out.sample_id,
        " -L chrM:8025-9144 "  
    )

    LIFTOVER_AND_COMBINE_VCFS(
        CALL_MUTECT_SHIFTED.out.vcf,
        CALL_MUTECT_STANDARD.out.vcf,
        ALIGN_STANDARD_MITO.out.sample_id,
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
        ALIGN_STANDARD_MITO.out.sample_id,
        blacklist,
        blacklist_index
    )

    SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES(
        mito_fasta,
        mito_index,
        mito_dict,
        INITIAL_FILTER.out.vcf,
        INITIAL_FILTER.out.tbi,
        ALIGN_STANDARD_MITO.out.sample_id
    )

    GET_CONTAMINATION(
        SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out.vcf
    )

    // CREATE_TABLE(GET_CONTAMINATION.out.minor_hg, GET_CONTAMINATION.out.major_hg)

    output(
        ALIGN_STANDARD_MITO.out.bam.concat(
            ALIGN_STANDARD_MITO.out.bai, 
            ALIGN_STANDARD_MITO.out.metrics,
            COLLECT_WGS_METRICS.out.sensitivity,
            COLLECT_WGS_METRICS.out.metrics,
            LIFTOVER_AND_COMBINE_VCFS.out.rejected,
            SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out.vcf,
            SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out.tbi,
            INITIAL_FILTER.out.vcf,
            INITIAL_FILTER.out.tbi,
            MERGE_STATS.out,
            GET_CONTAMINATION.out.major_hg,
            GET_CONTAMINATION.out.minor_hg,
            LIFTOVER_AND_COMBINE_VCFS.out.merged_vcf,
            // CREATE_TABLE.out,
            GET_CONTAMINATION.out.major_level,
            GET_CONTAMINATION.out.minor_level,
            GET_CONTAMINATION.out.contamination_file
        )
    )
}
