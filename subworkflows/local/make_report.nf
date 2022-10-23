//
// Aggregate all samples metrics in a single CSV file
//
include { AGGREGATE_SAMPLES_IN_CSV  } from '../../modules/local/custom/aggregate_samples_in_csv.nf'
include { CREATE_JSON               } from '../../modules/local/custom/create_sample_json.nf'

workflow make_report {
    take:
        contamination     // channel: [ val(meta), contam ]
        dup_metrics       // channel: [ val(meta), dup_metrics ]

    main:
        sample_metrics = contamination
            .join(dup_metrics)

        CREATE_JSON(sample_metrics)
        AGGREGATE_SAMPLES_IN_CSV(CREATE_JSON.out.collect())

    emit:
        csv = AGGREGATE_SAMPLES_IN_CSV.out  // channel: all_samples.csv
}
