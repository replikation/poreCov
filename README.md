![logo](data/logo/mobile_logo.png)
**poreCov | SARS-CoV-2 Workflow for nanopore sequencing data**   
===
![](https://img.shields.io/github/v/release/replikation/poreCov)
![](https://img.shields.io/badge/nextflow-20.10.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/uses-singularity-yellow.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4153510.svg)](https://doi.org/10.5281/zenodo.4153510)

![](https://github.com/replikation/nCov/workflows/Syntax_check/badge.svg)


[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 

* by Christian Brandt

> **Featured here:**
> Franziska Hufsky, Kevin Lamkiewicz, Alexandre Almeida, Abdel Aouacheria, Cecilia Arighi, Alex Bateman, Jan Baumbach, Niko Beerenwinkel, Christian Brandt, Marco Cacciabue, Sara Chuguransky, Oliver Drechsel, Robert D Finn, Adrian Fritz, Stephan Fuchs, Georges Hattab, Anne-Christin Hauschild, Dominik Heider, Marie Hoffmann, Martin Hölzer, Stefan Hoops, Lars Kaderali, Ioanna Kalvari, Max von Kleist, Renó Kmiecinski, Denise Kühnert, Gorka Lasso, Pieter Libin, Markus List, Hannah F Löchel, Maria J Martin, Roman Martin, Julian Matschinske, Alice C McHardy, Pedro Mendes, Jaina Mistry, Vincent Navratil, Eric P Nawrocki, Áine Niamh O’Toole, Nancy Ontiveros-Palacios, Anton I Petrov, Guillermo Rangel-Pineros, Nicole Redaschi, Susanne Reimering, Knut Reinert, Alejandro Reyes, Lorna Richardson, David L Robertson, Sepideh Sadegh, Joshua B Singer, Kristof Theys, Chris Upton, Marius Welzel, Lowri Williams, Manja Marz, Computational strategies to combat COVID-19: useful tools to accelerate SARS-CoV-2 and coronavirus research, Briefings in Bioinformatics, bbaa232, https://doi.org/10.1093/bib/bbaa232.

## What is this Repo?

* poreCov is a SARS-CoV-2 analysis workflow for nanopore data (via the ARTIC protocol) or just SARS-CoV-2 genomes (fasta) in general
    * ARTIC lab protocols and reference are [available here](https://artic.network/ncov-2019)
* QC, plots and overviews are included to decrease post analytic "downtime"
    * was the PCR coverage on each position enough?
    * is the quality of the genome good?
    * what variant mutations are present?
* is nanopore sequencing accurate enough for SARS-CoV-2 sequencing? [yes](https://www.nature.com/articles/s41467-020-20075-6)

## Quality Metrics (default)

* Regions with coverage of 20 or less are masked ("N")
* Genomequality is compared to NC_045512.2
    * `--rki` adds genome quality assessment based on [RKIBioinformaticsPipelines/president](https://gitlab.com/RKIBioinformaticsPipelines/president)
* Pangolin lineages are determined
* nextstrain clades are determined including mutation infos
* reads are classified to human and SARS-CoV-2 to check for possible contamination and sample prep issues

Table of Contents
=================

* [Installation](#Installation)
* [Run poreCov](#Run-poreCov)
    * [Example commands](#Example-commands)
    * [Help](#Help)
* [Workflow](#Workflow)
* [References and Metadata for tree construction](#References-and-Metadata-for-tree-construction)
    * [References](#References)
    * [Metadata](#Metadata)
* [Literature to cite](#Literature-to-cite)


# Installation

**Dependencies**

* one of these:
>   * docker
>   * singularity
>   * conda (install singularity + nextflow, usually not cluster compatible)

* these:
>   * nextflow + java runtime

* optional one of these if you want to start from basecalling (fast5)
>   * local guppy installation (see oxford nanopore installation guide)
>   * docker (with nvidia toolkit installed)
>   * singularity (with --nv support)
>   * dont have a gpu? all the above can be run with "cpu only" but takes ages

**Installation links**

* Docker installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce)
    * add docker to your User group via `sudo usermod -a -G docker $USER`
* Singularity installation [here](https://singularity.lbl.gov/install-linux)
    * if you cant use docker
* Conda installation [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
    * not natively integrated, you can do conda install singularity nextflow in a new environment and execute poreCov via `-profile local,singularity`
    * not cluster compatible
* Nextflow via: `curl -s https://get.nextflow.io | bash`
    * a java runtime is needed (e.g. `sudo apt install -y default-jre`)
    * move `nextflow` executable to a path location

# Run poreCov

* poreCov supports version control via `-r` this way you can run everything reproducable (e.g. `-r 0.7.2`)
* poreCov relases are listed [here](https://github.com/replikation/poreCov/releases)
* add `-r <version>` e.g. `nextflow run replikation/poreCov -r 0.7.0 -profile test_fastq,local,singularity`
* as nextflow manages the code you can update your "local poreCov" code via: `nextflow pull replikation/poreCov`

## Example commands

```bash
# Start from basecalling:
nextflow run replikation/poreCov --fast5 fast5/ -r 0.7.2 \
    --cores 32 -profile local,docker \
    --rki

# add multiple "combined" fastq.gz files (one sample per .fastq.gz file)
nextflow run replikation/poreCov -r 0.7.2 \
    --fastq 'samples/samplenumber_*.fastq.gz' \
    --cores 32  -profile local,docker

# use the fastq_pass dir from "guppy/basecalling output"
# this dir schould contain "barcode??/" dirs   e.g. fastq_pass/barcode01/ fastq_pass/barcode02/
nextflow run replikation/poreCov --fastq_pass 'fastq_pass/' -r 0.7.2 \
    --cores 32  -profile local,docker

# rename barcodes automatically by providing an input file, also using another primer scheme
nextflow run replikation/poreCov --fast5 fast5_dir/ --samples sample_names.csv \
    --extended --primerV V1200 --output results --rki -profile local,docker
```

## Help

* workflows and inputs are described here:

```bash
nextflow run replikation/poreCov -r 0.7.2 --help
# or git clone and
./poreCov.nf --help
```

# Workflow

* poreCov was coded with "easy to use" in mind, while staying flexible
* therefor we provide a few input types which adjusts the workflow automatically (see image below)
  * fast5 raw data, fastq files (one sample per file), fastq_pass (the basecalling output) or fasta (supports multifastas)
* primer schemes for ARTIC can be V1, V2, V3(default) or V1200 (the 1200bp amplicon ones)

![workflow](data/figures/workflow.png)


# Literature / References to cite
For citing etc. check out these programs used for poreCov:
* [nextflow](https://www.nextflow.io/index.html)
* [artic protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)
* [pangolin](https://github.com/hCoV-2019/pangolin)
* [medaka](https://github.com/nanoporetech/medaka)
* [president](https://gitlab.com/RKIBioinformaticsPipelines/president)
* [nextclade](https://clades.nextstrain.org/)
* [kraken2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)
* [krona](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385)
