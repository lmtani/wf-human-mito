# Mitochondria Workflow

This repository contains a [nextflow](https://www.nextflow.io/) workflow
for running mitochondrial analysis. This workflows is heavily inspired by the [gatk-workflows/gatk4-mitochondria-pipeline](https://github.com/gatk-workflows/gatk4-mitochondria-pipeline) and by
the [epi2me-labs/wf-template](https://github.com/epi2me-labs/wf-template).


## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or conda is installed.


**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run lmtani/wf-human-mito -r main --help
```

to see the options for the workflow.

**Workflow outputs**

WIP


## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
