![logo](data/logo/mobile_logo.png)
**nCov19 Workflow for nanopore sequencing data**   
===

![](https://img.shields.io/badge/nextflow-20.01.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/uses-singularity-yellow.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)


![](https://github.com/replikation/nCov/workflows/Syntax_check/badge.svg)


[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 

* by Christian Brandt
* **under development** 

## What is this Repo?

* a general nCov analysis workflows using nanopore data with the ARTIC and Augur tools
    * ARTIC protocol and reference from [here](https://artic.network/ncov-2019)
* a few handy QC and plots are included to decrease post analytic "downtime"
    * e.g. was the PCR coverage on each position enough?

## Workflow

* nCov was coded with "easy to use" in mind, while staying flexible
* the default use case is mostly fast5 raw-data to "results"
* however by providing fastq or fasta instead the workflow skips over the corresponding modules (see workflow figure below)
* only by providing `--metadata` tree construction starts

![workflow](data/figures/workflow.png)


## References and Metadata for tree construction
### References
* by default the nCov Workflow uses a few hundred nCov strains automatically to build the tree
    * these files are from ENA
* this behaviour can be "replaced" by providing a multifasta reference file via `--references`
* please make sure that your fasta file header does not contain "strange" symbols like `" '  | / \ :`
    * this causes issues in some of the tools used here

### Metadata
* for tree constructions metad ata is mandatory (time and location)
    * location can be set to "unknown"
* provide a metadata file via `--metadata` to create a correct tree
    * if you use you `--references` add this information to you metadata file too
* style of the meta data (tab separated):
    * strain is the name of the fasta header without `>`
    * date has to be in the format `YYYY-MM-DD` or `YYYY-MM`

```csv
strain	country	date
LC528232-unknown-2020-02-10	unknown	2020-02-10
LC528233-unknown-2020-02-10	unknown	2020-02-10
LC534418-Japan-2020-02-14	Japan	2020-02-14
LC534419-Japan-2020-03-09	Japan	2020-03-09
MN908947-China-Dec-2019	China	2019-12
```

## Installation

**Dependencies**

* one of these:
>   * docker
>   * singularity
>   * conda (NOT YET IMPLEMENTED)
* all of these:
>   * a local guppy installation if you dont have a gpu docker
>      * not needed if you only use fastq or fasta as input
>   * nextflow + java runtime 

* Docker installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce)
    * add docker to your User group via `sudo usermod -a -G docker $USER`
* Singularity installation [here](https://singularity.lbl.gov/install-linux)
* Conda installation [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
* Nextflow installation [here](https://www.nextflow.io/)
* move or add the nextflow executable to a bin path

## Current implementations and a few planned features

* [x] basecalling via guppy (either local install or gpu docker)
    * [x] single sample input (fast5)
    * [x] barcoded fast5 input
* [x] [artic](https://github.com/artic-network/fieldbioinformatics) via nanopore
    * [x] read coverage plots to validate PCR product coverage for genome construction
    * [ ] quality of genomes summary
* [x] Reference based sample comparision
    * [x] augur based, see also [Link](https://nextstrain.org/help/coronavirus/SARS-CoV-2)
    * [x] timetree / location calc. inclusion
* [x] toytree vis plot
* [x] docker engine
* [x] singularity engine
* [ ] conda

## Help

* workflows and inputs are described here:

```bash
nextflow run replikation/nCov --help
# or
./nCov.nf --help
```
