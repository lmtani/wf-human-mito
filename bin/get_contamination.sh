#!/usr/bin/env bash

set -e

input_vcf=$1


PARENT_DIR="$(dirname "${input_vcf}")"
java -jar /usr/mtdnaserver/haplocheckCLI.jar "${PARENT_DIR}"

sed 's/\"//g' output > output-noquotes
