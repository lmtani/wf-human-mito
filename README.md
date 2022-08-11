# Human Mitochondrial Workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow for running mitochondrial analysis. This workflow is heavily inspired by the [gatk-workflows/gatk4-mitochondria-pipeline](https://github.com/gatk-workflows/gatk4-mitochondria-pipeline) and by the [nf-core](https://nf-co.re) community.

## 🔧 Setup

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or [conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of the required software. Both methods are automated out-of-the-box provided either docker or conda is installed.

<p align="center">
    <img title="Workflow diagram" src="./docs/wf-human-mito-metro.png" width=75%>
</p>

## 📥 Inputs

- Pairs of FASTQ file. One pair for each sample or
- Alignment file, one per sample. Accepts BAM or CRAM formats.
- Human Genome Reference - choose the one that best suits your needs.

### Notes

- Allows paired FASTQs, alignments or both.
- Some of the human genome can be downloaded [here](https://cloud.google.com/life-sciences/docs/resources/public-datasets/reference-genomes), hosted by Google. Example: [GRCh38](https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0). This workflow needs to have access to *.fasta, .dict, .fai, 64.ann, 64.amb, 64.sa, 64.pac, 64.alt* files.
- The *--reference* parameter must point to the same reference of the input alignments.

## ⚙ Running

```bash
# For help:
nextflow run lmtani/wf-human-mito -r main --help

# Example:
nextflow run lmtani/wf-human-mito -r main \
    --fastq 'fastqs/*_R{1,2}.fq.gz' \ # ex: sample_R1.fq.gz and sample_R2.fq.gz
    --alignments 'bams/*.cra{m,i}' \  # ex: sample.cram and sample.crai
    --reference /refs/Homo_sapiens_assembly38.fasta \
    --outdir outdir \
    -profile conda   # or docker
```

## 📤 Outputs

- Alignment in BAM format (*outdir/alignments/*)
- Variants in VCF format (*outdir/variants/*)
- CSV file with informations, e.g: Haplotype groups (major and minor), coverage, etc. for all samples.
- All intermediate outputs (ex: *outdir/workspace/align_raw_reads/* contains all whole genome alignments)

## Useful links

- [nextflow](https://www.nextflow.io/)
- [docker](https://www.docker.com/products/docker-desktop)
- [conda](https://docs.conda.io/en/latest/miniconda.html)

### TODOs

- [ ] In `COLLECT_WGS_METRICS` process, the input value "readLen" should be dynamic. This info is inside COLLECT_ALIGNMENT_METRICS output.
