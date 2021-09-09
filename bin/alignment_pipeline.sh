#!/usr/bin/env bash

set -ex

reference=$1
fq1=$2
fq2=$3
basename=$4


id=$(zcat "$fq1" | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g' | cut -d "-" -f 2)
platform_unit=$(zcat "$fq1" | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g' | cut -d "-" -f 1)
sm=$(zcat "$fq1" | head -n 1 | grep -Eo "[ATGCN+]+$")


bwa mem -K 100000000 -v 3 -t 4 -Y \
    -R $(echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPU:$platform_unit\tPL:ILLUMINA") \
    "$reference" \
    "$fq1" "$fq2" | samtools sort -o algn.bam -
samtools index algn.bam


picard \
    MarkDuplicates \
    INPUT=algn.bam \
    OUTPUT="$basename.markdup.bam" \
    METRICS_FILE="$basename.dup.metrics" \
    VALIDATION_STRINGENCY=SILENT \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    ASSUME_SORT_ORDER="queryname" \
    CLEAR_DT="false" \
    ADD_PG_TAG_TO_READS=false

picard \
    SortSam \
    INPUT="$basename.markdup.bam" \
    OUTPUT="$basename.sorted.bam" \
    SORT_ORDER="coordinate" \
    CREATE_INDEX=true \
    MAX_RECORDS_IN_RAM=300000
