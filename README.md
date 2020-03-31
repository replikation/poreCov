![logo](data/logo/mobile_logo.png)
**nCov Workflows**   |   Institut für Infektionsmedizin und Krankenhaushygiene
===

![](https://img.shields.io/badge/nextflow-20.01.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)
![](https://github.com/replikation/nCov/workflows/Syntax_check/badge.svg)


[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 

* by Christian Brandt
* **under heavy development**

## What is this Repo?

* general nCov analysis workflows collection for the JUH  

## Workflows

+ for live basecalling use [this repos](https://github.com/replikation/docker_pipelines)
+ currently: 
    + use artic workflow for genome construction 
    + use augur to analyse nCov19 genomes afterwards (either from artic results or use you own fastas)

### Current implementations and a few planned features

* [x] [artic](https://github.com/artic-network/fieldbioinformatics) via nanopore
    * [x] single flongle input (one fastq file of a sample)
    * [ ] demultiplex, barcoded fastq input
* [ ] metagenomic approach via nanopore
* [ ] direct DNA sequencing via nanopore
* [ ] direct RNA sequencing via nanopore
* [x] Reference based sample comparision
    * [x] mafft based comparison (WIP time interference is missing here - use augur for now)
    * [x] augur based comparision, see also [Link](https://nextstrain.org/help/coronavirus/SARS-CoV-2)
    * [x] time / location calc. inclusion
* [x] toytree vis


## Installation

**Dependencies**

>   * docker (add docker to your Usergroup, so no sudo is needed)
>   * nextflow + java runtime 
>   * git (should be already installed)
>   * wget (should be already installed)

* Docker installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce)
* Nextflow installation [here](https://www.nextflow.io/)
* move or add the nextflow executable to a bin path
* add docker to your User group via `sudo usermod -a -G docker $USER`

## Help

* workflows and inputs are described here:

```bash
nextflow run replikation/nCov --help
# or
./nCov.nf --help
```