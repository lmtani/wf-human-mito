#!/usr/bin/env bash

set -ex

f1=$1
f2=$2
sample_id=$3


READ_NAME=$(zcat "$f1" | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g' | cut -d "-" -f 1)
RUN_ID=$(echo "$READ_NAME" | awk -F "_" '{print $1"_"$2"_"$3"_"$4}')
PLATFORM_UNIT=$(echo "$READ_NAME" | awk -F "_" '{print $2}')

picard FastqToSam \
    FASTQ="$f1" \
    FASTQ2="$f2" \
    OUTPUT="$sample_id.unmaped.bam" \
    READ_GROUP_NAME=A \
    SAMPLE_NAME="$sample_id" \
    LIBRARY_NAME="$RUN_ID" \
    PLATFORM_UNIT="$PLATFORM_UNIT" \
    PLATFORM=illumina \
    SEQUENCING_CENTER=BI
