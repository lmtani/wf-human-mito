---
layout: default
title: Quickstart
nav_order: 1
---

# 🚀 Quickstart

This is a workflow (written in Nextflow) to analyse mitochondrial sequencing experiments. Is is heavily inspired by the [gatk-workflows/gatk4-mitochondria-pipeline](https://github.com/gatk-workflows/gatk4-mitochondria-pipeline) and by the [epi2me-labs/wf-template](https://github.com/epi2me-labs/wf-template).


1. Install Nextflow in your system (see [official documentation](https://www.nextflow.io/docs/latest/getstarted.html#installation))
2. Download the GRCh38 release of the Human Genome and its index files. You need to save all these files in the same directory.

    ??? summary "Links from [Google Cloud Life Sciences](https://cloud.google.com/life-sciences/docs/resources/public-datasets/reference-genomes)"
        - [Homo_sapiens_assembly38.fasta](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta)
        - [Homo_sapiens_assembly38.dict](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dict)
        - [Homo_sapiens_assembly38.fasta.fai ](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai)
        - [Homo_sapiens_assembly38.fasta.64.alt](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt)
        - [Homo_sapiens_assembly38.fasta.64.amb](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb)
        - [Homo_sapiens_assembly38.fasta.64.ann](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann)
        - [Homo_sapiens_assembly38.fasta.64.bwt](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt)
        - [Homo_sapiens_assembly38.fasta.64.pac](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac)
        - [Homo_sapiens_assembly38.fasta.64.sa](https://storage.cloud.google.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa)

3. Put your FASTQ files (or symbolic links `ln -s <original> <link>`) in one directory. Use a name pattern in a way that you can identify your samples from the prefix.

    ??? summary "Example"
        ```bash
        $ ls test_data/*.fq.gz
        test_data/SAMPLE-A_R1.fq.gz  test_data/SAMPLE-C_R1.fq.gz
        test_data/SAMPLE-A_R2.fq.gz  test_data/SAMPLE-C_R2.fq.gz
        test_data/SAMPLE-B_R1.fq.gz  test_data/SAMPLE-D_R1.fq.gz
        test_data/SAMPLE-B_R2.fq.gz  test_data/SAMPLE-D_R2.fq.gz
        ```

4. Run the Human Mitochondrial Workflow

    ??? summary "Example"
        ```bash
        ./nextflow run lmtani/wf-human-mito \
            --fastq "test_data/*_R{1,2}.fq.gz" \
            --reference /path/to/Homo_sapiens_assembly38.fasta \
            --outdir name-your-output-directory
        ```

In the end, you will have variants (VCF) and alignment (BAM) for each sample and a CSV file with information about each sample (coverage, haplogroup, etc). For more details about the meaning of each column, please see the [📦 Outputs section](outputs)

??? summary "Example"
    ```bash
    $ ls -1 name-your-output-directory/
    alignments/            # Alignments (BAM)
    all_samples.csv        # Information about each sample
    execution/             # Pipeline execution details
    variants/              # Variants (VCF)
    workspace/             # Outputs of each process in the pipeline
    ```
