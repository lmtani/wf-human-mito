//
// Call mitochondrial variants, using shifted approach
// described in https://gnomad.broadinstitute.org/news/2020-11-gnomad-v3-1-mitochondrial-dna-variants/#mtdna-calling-pipeline-for-single-samples
//
include {
          CALL_VARIANTS as CALL_DEFAULT;
          CALL_VARIANTS as CALL_SHIFTED        } from '../../subworkflows/local/mutect2_variant_call.nf'
include { GATK4_FILTERMUTECTCALLS              } from '../../modules/local/gatk4/filtermutectcalls/main'
include { GATK4_LEFTALIGNANDTRIMVARIANTS       } from '../../modules/nf-core/gatk4/leftalignandtrimvariants/main'
include { GATK4_MERGEMUTECTSTATS               } from '../../modules/nf-core/gatk4/mergemutectstats/main'
include { GATK4_SELECTVARIANTS                 } from '../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_VARIANTFILTRATION              } from '../../modules/nf-core/gatk4/variantfiltration/main'
include { HAPLOCHECK                           } from '../../modules/nf-core/haplocheck/main'
include { MERGE_VCFS                           } from '../../modules/local/picard/merge_vcfs.nf'
include { PICARD_LIFTOVERVCF                   } from '../../modules/nf-core/picard/liftovervcf/main'
include { MOSDEPTH                             } from '../../modules/nf-core/mosdepth/main'

workflow variant_call {
    take:
        reads           // channel: [ val(meta), ubam ]
        standard_genome // channel: [ fasta, dict, index, amb, ann, bwt, pac, sa, alt, intervals ]
        shifted_genome  // channel: [ fasta, dict, index, amb, ann, bwt, pac, sa, alt, intervals ]
    main:

        // To gather all QC reports for MultiQC
        ch_reports  = Channel.empty()

        // To gather used softwares versions for MultiQC
        ch_versions = Channel.empty()

        CALL_DEFAULT(reads, " -L chrM:576-16024 ", "standard", standard_genome)
        ch_versions = ch_versions.mix(CALL_DEFAULT.out.versions)

        CALL_SHIFTED(reads, " -L chrM:8025-9144 ", "shifted", shifted_genome)
        ch_versions = ch_versions.mix(CALL_SHIFTED.out.versions)

        PICARD_LIFTOVERVCF(CALL_SHIFTED.out.vcf, params.genome.mito_dict, params.genome.shift_back_chain, params.genome.mito_fasta)
        ch_versions = ch_versions.mix(PICARD_LIFTOVERVCF.out.versions)

        MOSDEPTH(CALL_DEFAULT.out.alignment.map { it -> [ it[0], it[1], it[2] ]}, [], params.genome.mito_fasta)
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
        mosdepth_metrics = MOSDEPTH.out.global_txt
                                       .join(MOSDEPTH.out.summary_txt)
                                       .join(MOSDEPTH.out.per_base_bed)
                                       .join(MOSDEPTH.out.per_base_csi)
        ch_reports = ch_reports.mix(mosdepth_metrics.map { it -> [ it[1], it[2], it[3], it[4] ]})

        vcfs_channel = CALL_DEFAULT.out.vcf.join(PICARD_LIFTOVERVCF.out.vcf_lifted)
        MERGE_VCFS(vcfs_channel)
        ch_versions = ch_versions.mix(MERGE_VCFS.out.versions)

        stats_channel = CALL_DEFAULT.out.mutect_stats.join(CALL_SHIFTED.out.mutect_stats).map { it -> [ it[0], [ it[1], it[2] ] ]}
        GATK4_MERGEMUTECTSTATS(stats_channel)
        ch_versions = ch_versions.mix(GATK4_MERGEMUTECTSTATS.out.versions)

        filter_variants_channel = MERGE_VCFS.out.vcf.join(MERGE_VCFS.out.idx).join(GATK4_MERGEMUTECTSTATS.out.stats).map( it ->
            [ it[0], it[1], it[2], it[3], [], [], [], [] ]
        )
        GATK4_FILTERMUTECTCALLS(filter_variants_channel, params.genome.mito_fasta, params.genome.mito_index, params.genome.mito_dict)
        ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions)

        variantfiltration_channel = GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi)
        GATK4_VARIANTFILTRATION(variantfiltration_channel, params.genome.mito_fasta, params.genome.mito_index, params.genome.mito_dict)
        ch_versions = ch_versions.mix(GATK4_VARIANTFILTRATION.out.versions)

        left_align_channel = GATK4_VARIANTFILTRATION.out.vcf.join(GATK4_VARIANTFILTRATION.out.tbi).map( it -> [ it[0], it[1], it[2], []] )
        GATK4_LEFTALIGNANDTRIMVARIANTS(left_align_channel, params.genome.mito_fasta, params.genome.mito_index, params.genome.mito_dict)
        ch_versions = ch_versions.mix(GATK4_LEFTALIGNANDTRIMVARIANTS.out.versions)

        GATK4_SELECTVARIANTS(GATK4_LEFTALIGNANDTRIMVARIANTS.out.vcf.join(GATK4_LEFTALIGNANDTRIMVARIANTS.out.tbi))
        ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions)

        HAPLOCHECK(GATK4_SELECTVARIANTS.out.vcf)
        ch_versions = ch_versions.mix(HAPLOCHECK.out.versions)

        ch_reports = ch_reports
                            .mix(CALL_DEFAULT.out.dup_metrics.map { it -> [ it[1] ]})
                            .mix(CALL_DEFAULT.out.txt_metrics.map{ it -> it[1] })
                            .mix(HAPLOCHECK.out.txt.map{ it -> it[1] })

    emit:
        alignment          = CALL_DEFAULT.out.alignment                                                // channel: [ val(sample_id), bam, bai ]
        contamination      = HAPLOCHECK.out.txt                                                        // channel: [ val(sample_id), contam ]
        mutect_vcf         = GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi)     // channel: [ val(sample_id), vcf, tbi ]
        reports            = ch_reports
        versions           = ch_versions
        multiqc_rename     = CALL_DEFAULT.out.multiqc_rename
}
