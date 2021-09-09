#!/usr/bin/env bash

set -e

f1=$1
f2=$2
sample_id=$3


SAMPLE_NAME=$(zcat "$f1" | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g' | cut -d "-" -f 2)
PLATFORM_UNIT=$(zcat "$f1" | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g' | cut -d "-" -f 1)
BARCODE=$(zcat "$f1" | head -n 1 | grep -Eo "[ATGCN+]+$")


picard FastqToSam \
    FASTQ="$f1" \
    FASTQ2="$f2" \
    OUTPUT="$sample_id.unmaped.bam" \
    READ_GROUP_NAME=A \
    SAMPLE_NAME="$SAMPLE_NAME" \
    LIBRARY_NAME="$SAMPLE_NAME-$BARCODE" \
    PLATFORM_UNIT="$PLATFORM_UNIT" \
    PLATFORM=illumina \
    SEQUENCING_CENTER=BI
