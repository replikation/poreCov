![logo](data/logo/mobile_logo.png)
**poreCov | SARS-CoV-2 Workflow for nanopore sequencing data**   
===

![](https://img.shields.io/badge/nextflow-20.07.0-brightgreen)
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

* poreCov is a general nCov analysis workflows using nanopore data with the ARTIC and Augur tools
    * ARTIC protocol and reference from [here](https://artic.network/ncov-2019)
* a few handy QC and plots are included to decrease post analytic "downtime"
    * e.g. was the PCR coverage on each position enough?
* is nanopore sequencing accurate enough for SARS-CoV-2 sequencing? [yes](https://www.nature.com/articles/s41467-020-20075-6)

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
>   * conda (NOT YET IMPLEMENTED)

* these:
>   * a local guppy installation or a docker with gpu support
>      * not needed if you use fastq or fasta as input
>   * nextflow + java runtime 

**Installation links**

* Docker installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce)
    * add docker to your User group via `sudo usermod -a -G docker $USER`
* Singularity installation [here](https://singularity.lbl.gov/install-linux)
    * if you cant use docker
* Conda installation [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
    * not natively integrated, you can do conda install singularity nextflow in a new environment and execute poreCov via `-profile local,singularity`
* Nextflow via: `curl -s https://get.nextflow.io | bash`
    * a java runtime is needed
    * move `nextflow` executable to a path location

# Run poreCov

## Example commands

```bash
# just do basecalling and assembly with QC / lineage:
nextflow run replikation/poreCov --dir fast5/ \
    --cores 32 -profile local,docker

# use fastq files with workflow provided references for tree construction
nextflow run replikation/poreCov  \
    --fastq 'samples/samplenumber_*.fastq.gz' --metadata metadata.csv \
    --cores 32  -profile local,docker

# build a tree by using your own references
nextflow run replikation/poreCov --fasta 'nCov-genomes/*.fasta' \
    --metadata metadata.csv --references references_multifasta.fa \
    --cores 32  -profile local,docker
```

## Help

* workflows and inputs are described here:

```bash
nextflow run replikation/poreCov --help
# or git clone and
./poreCov.nf --help
```

# Workflow

* poreCov was coded with "easy to use" in mind, while staying flexible
* the default use case is fast5 raw-data to "results"
* however by providing fastq or fasta instead poreCov skips over the corresponding steps
* provide `--metadata` to start the tree construction
* genomes with 3 or more `N's` are excluded (not ignoring the 70bp at N and C-terminal)
  * can be modified via `--rm_N_genome`
* primer schemes for ARTIC can be V1,V2,V3 or the 1200bp ones (see help)

![workflow](data/figures/workflow.png)

# References and Metadata for tree construction
## References
* by default the nCov Workflow uses a few hundred nCov strains automatically to build the tree
    * these files are from ENA
    * metadata for this is automatically added
* this behaviour can be "replaced" by providing a multifasta reference file via `--references`
* please make sure that these fastaheader dont contain "strange" symbols like `" '  | / \ : & $`
    * this causes issues in some of the tools used here

## Metadata
* for tree constructions `--metadata file.tsv` is mandatory
* style of the metadata: 
    * tab separated: `\t`
    * header: `strain   time    location`
    * one sample per line
    * **strain**: 
        * name of the fastq file without any suffix (`--fastq`)
        * barcode01 barcode02 etc. for fast5 data inside a directory (`--dir`)
        * fasta header name (without >) for one genome (`--fasta`) or multiple references (`--reference`)
    * **location** can be set to "unknown"
    * **date** format is `YYYY-MM-DD` or `YYYY-MM` or `YYYY`
* example:

```csv
strain	country	date
LC528232-unknown-2020-02-10	unknown	2020-02-10
LC528233-unknown-2020-02-10	unknown	2020-02-10
LC534418-Japan-2020-02-14	Japan	2020-02-14
LC534419-Japan-2020-03-09	Japan	2020-03-09
MN908947-China-Dec-2019	China	2019-12
barcode01   Germany 2020-01-01
```

# Literature / References to cite
For citing etc. check out these programs used for poreCov:
* [nextflow](https://www.nextflow.io/index.html)
* [artic protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)
* [augur](https://github.com/nextstrain/augur)
* [pangolin](https://github.com/hCoV-2019/pangolin)
    * [iqtree](http://www.iqtree.org/#download)
    * [mafft](https://mafft.cbrc.jp/alignment/software/)
    * [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)
* [toytree](https://github.com/eaton-lab/toytree)
* [medaka](https://github.com/nanoporetech/medaka)
* [president](https://gitlab.com/RKIBioinformaticsPipelines/president)
