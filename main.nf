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
    BWA_ALIGN as ALIGN_SHIFTED_MITO;
    BWA_ALIGN as ALIGN_STANDARD_MITO;
    BWA_ALIGN_FROM_UBAM;
    CALL_MUTECT as CALL_MUTECT_SHIFTED;
    CALL_MUTECT as CALL_MUTECT_STANDARD;
    COLLECT_ALIGNMENT_METRICS;
    COLLECT_WGS_METRICS;
    CREATE_ALL_SAMPLES_CSV;
    CREATE_JSON;
    FASTQ_TO_UBAM;
    FILTER as INITIAL_FILTER;
    GET_CONTAMINATION;
    LIFTOVER_AND_COMBINE_VCFS;
    MERGE_STATS
    SELECT_MITO_READS;
    SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES} from './modules.nf'

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


workflow separate_mitochondrion {
    take: reads
    main:

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

        FASTQ_TO_UBAM(reads)
        sample_id = FASTQ_TO_UBAM.out.sample_id
        BWA_ALIGN_FROM_UBAM(
            sample_id,
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
            sample_id,
            BWA_ALIGN_FROM_UBAM.out.bam,
            BWA_ALIGN_FROM_UBAM.out.bai,
            human_fasta,
            human_dict,
            human_index,
        )
    emit:
        ubam = SELECT_MITO_READS.out.ubam
        sample_id = sample_id
}


workflow variant_call {
    take: 
        ubam
        sample_id
    main:
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
    
        ALIGN_STANDARD_MITO(
                ubam,
                sample_id,
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
            ubam,
            sample_id,
            shifted_fasta,
            shifted_dict,
            shifted_index,
            shifted_amb,
            shifted_ann,
            shifted_bwt,
            shifted_pac,
            shifted_sa
        )


        alignment = ALIGN_STANDARD_MITO.out.bam
        alignment_index = ALIGN_STANDARD_MITO.out.bai
        COLLECT_ALIGNMENT_METRICS(
            alignment,
            alignment_index,
            mito_fasta,
            mito_dict,
            mito_index
        )

        COLLECT_WGS_METRICS(
            alignment,
            alignment_index,
            mito_fasta,
            sample_id,
            300
        )

        CALL_MUTECT_STANDARD(
            alignment,
            alignment_index,
            mito_fasta,
            mito_dict,
            mito_index,
            "standard",
            sample_id,
            " -L chrM:576-16024 "
        )

        CALL_MUTECT_SHIFTED(
            ALIGN_SHIFTED_MITO.out.bam,
            ALIGN_SHIFTED_MITO.out.bai,
            shifted_fasta,
            shifted_dict,
            shifted_index,
            "shifted",
            sample_id,
            " -L chrM:8025-9144 "  
        )

        LIFTOVER_AND_COMBINE_VCFS(
            CALL_MUTECT_SHIFTED.out.vcf,
            CALL_MUTECT_STANDARD.out.vcf,
            sample_id,
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
            sample_id,
            blacklist,
            blacklist_index
        )

        SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES(
            mito_fasta,
            mito_index,
            mito_dict,
            INITIAL_FILTER.out.vcf,
            INITIAL_FILTER.out.tbi,
            sample_id,
        )

        GET_CONTAMINATION(
            SPLIT_MULTIALLELICS_AND_REMOVE_NON_PASS_SITES.out.vcf
        )
    emit:
        dup_metrics = ALIGN_STANDARD_MITO.out.metrics
        wgs_metrics = COLLECT_WGS_METRICS.out.metrics
        contam_metrics = GET_CONTAMINATION.out.contamination_file
        algn_metrics = COLLECT_ALIGNMENT_METRICS.out
        bam = ALIGN_STANDARD_MITO.out.bam
        bai = ALIGN_STANDARD_MITO.out.bai
        vcf = INITIAL_FILTER.out.vcf
        tbi = INITIAL_FILTER.out.tbi
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

    if (!params.reference) {
        helpMessage()
        println("")
        println("`--reference` is required")
        exit 1
    }


    reads = Channel.fromFilePairs("${params.fastq}", glob: true)

    separate_mitochondrion(reads)

    variant_call(separate_mitochondrion.out.ubam, separate_mitochondrion.out.sample_id)
  
    CREATE_JSON(
        variant_call.out.dup_metrics,
        variant_call.out.wgs_metrics,
        variant_call.out.contam_metrics,
        variant_call.out.algn_metrics,
        separate_mitochondrion.out.sample_id
    )

    CREATE_ALL_SAMPLES_CSV(
        CREATE_JSON.out.collect()
    )

    output(
        variant_call.out.bam.concat(
            variant_call.out.bai, 
            variant_call.out.vcf,
            variant_call.out.tbi,
            CREATE_ALL_SAMPLES_CSV.out
        )
    )
}
