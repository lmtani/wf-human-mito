# Human Mitochondrial Workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow for running mitochondrial analysis. This workflow is heavily inspired by the [gatk-workflows/gatk4-mitochondria-pipeline](https://github.com/gatk-workflows/gatk4-mitochondria-pipeline) and by the [nf-core](https://nf-co.re) community.

## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or conda is installed.

## Inputs

- Pairs of FASTQ file. One pair for each sample.
- Reference for human genome (GRCh38). [Files are available here](https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0).
  - .fasta, .dict, .fai, .ann, .amb, .sa, .pac, .alt

## Workflow options

To obtain the workflow, having installed `nextflow`, users can run:

```bash
nextflow run lmtani/wf-human-mito -r main --help
```

to see the options for the workflow.

## Workflow outputs

- Alignment in BAM format
- Variants in VCF format
- JSON with informations, e.g: Haplotype groups (major and minor), coverage, etc.

## Useful links

- [nextflow](https://www.nextflow.io/)
- [docker](https://www.docker.com/products/docker-desktop)
- [conda](https://docs.conda.io/en/latest/miniconda.html)

### TODOs

- [ ] In `COLLECT_WGS_METRICS` process, the input value "readLen" should be dynamic. This info is inside COLLECT_ALIGNMENT_METRICS output.
