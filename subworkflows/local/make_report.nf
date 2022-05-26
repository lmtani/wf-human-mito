//
// Aggregate all samples metrics in a single CSV file
//
include { CREATE_JSON               } from '../../modules/local/custom/create_sample_json.nf'
include { AGGREGATE_SAMPLES_IN_CSV  } from '../../modules/local/custom/aggregate_samples_in_csv.nf'

workflow make_report {
    take:
        contamination     // channel: [ val(sample_id), contam ]
        alignment_metrics // channel: [ val(sample_id), algn_metrics, theoretical_sensitivity ]
        alignment_wgs     // channel: [ val(sample_id), wgs_metrics ]
        dup_metrics       // channel: [ val(sample_id), dup_metrics ]

    main:
        sample_metrics = contamination
            .join(alignment_metrics)
            .join(alignment_wgs)
            .join(dup_metrics)

        CREATE_JSON(sample_metrics)
        AGGREGATE_SAMPLES_IN_CSV(CREATE_JSON.out.collect())

    emit:
        csv = AGGREGATE_SAMPLES_IN_CSV.out  // channel: all_samples.csv
}
