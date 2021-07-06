<p align="center">
  <img src="data/logo/mobile_logo.png" width="800" title="Workflow">
</p>

**poreCov | SARS-CoV-2 Workflow for nanopore sequencing data**   
===
![](https://img.shields.io/github/v/release/replikation/poreCov)
![](https://img.shields.io/badge/nextflow-20.10.0-brightgreen)
![](https://img.shields.io/badge/uses-Docker-blue.svg)
![](https://img.shields.io/badge/uses-Singularity-yellow.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)
![](https://github.com/replikation/poreCov/actions/workflows/nextflow-test.yml/badge.svg)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4153509.svg)](https://doi.org/10.5281/zenodo.4153509)

[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 


**Citation:**  
> **poreCov - an easy to use, fast, and robust workflow for SARS-CoV-2 genome reconstruction via nanopore sequencing**  
> Christian Brandt, Sebastian Krautwurst, Riccardo Spott, Mara Lohde, Mateusz Jundzill, Mike Marquet, Martin Hölzer  
> https://www.biorxiv.org/content/10.1101/2021.05.07.443089v1.

## What is this Repo?

* poreCov is a SARS-CoV-2 analysis workflow for nanopore data (via the [ARTIC protocol](https://artic.network/ncov-2019)) or SARS-CoV-2 genomes (fasta)
* the workflow is pre configured to simplify [data analysis](https://htmlpreview.github.io/?https://github.com/replikation/poreCov/blob/master/data/figures/index.html): 
<p align="left">
    <a href="https://htmlpreview.github.io/?https://github.com/replikation/poreCov/blob/master/data/figures/index.html">
        <img src="data/figures/report_summary.png" width="500" title="Report file">
</p>

Table of Contents
=================
<!--ts-->
* [1. Quick Setup (Ubuntu)](#1-quick-setup-ubuntu)
    * [1.1 Nextflow](#11-nextflow-the-workflow-manager)
    * [1.2 Container](#12-container-choose-one---they-manage-all-the-tools)
    * [1.3 Basecalling (optional)](#13-basecalling-optional)
* [2. Run poreCov](#2-run-porecov)
    * [2.1 Test run](#21-test-run)
    * [2.2 Quick run examples](#22-quick-run-examples)
    * [2.3 Extended Usage](#23-extended-usage)
* [3. Quality Metrics (default)](#3-quality-metrics-default)
* [4. Workflow](#4-workflow)
* [5. Literature / References to cite](#5-literature--references-to-cite)
* [6. Troubleshooting](#6-troubleshooting)
* [7. Time to results](#7-time-to-results)
<!--te-->

  # 1. Quick Setup (Ubuntu)
## 1.1 Nextflow (the workflow manager)
* poreCov needs [Nextflow](https://www.nextflow.io/index.html) and java run time (default-jre)
    * install java run time via:  `sudo apt install -y default-jre`
    * install Nextflow via:  `curl -s https://get.nextflow.io | bash && sudo mv nextflow /bin && sudo chmod 770 /bin/nextflow`
## 1.2 Container (choose one - they manage all the tools)
### Docker
* installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce) (recommended), alternatively via: `sudo apt install -y docker`
* add Docker to the user: `sudo usermod -a -G docker $USER`
### Singularity
* Singularity installation [here](https://singularity.lbl.gov/install-linux)
* if you can't use Docker
### Conda (not recommended)
* Conda installation [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
* install Nextflow and Singularity via conda (not cluster compatible) - and use the `singularity` profile
## 1.3 Basecalling (optional)
* only important if you want to do basecalling via GPU with the workflow:
    * local guppy installation (see oxford nanopore installation guide)
    * or: install nvidia Docker tool kit
    * or: Singularity (with --nv support)


# 2. Run poreCov
## 2.1 Test run
* validate your installation via test data:

```bash
# for a Docker installation
nextflow run replikation/poreCov -profile test_fastq,local,docker -r 0.8.0

# or for Singularity or conda installation
nextflow run replikation/poreCov -profile test_fastq,local,singularity -r 0.8.0
```

## 2.2 Quick run examples

* poreCov with basecalling and Docker
    * `--update` tryies to force the most recent lineage release version (optional)
    * `-r 0.8.0` specifies the workflow release from [here](https://github.com/replikation/poreCov/releases)
```bash
nextflow run replikation/poreCov --fast5 fast5/ -r 0.8.0 \
    --cores 6 -profile local,docker --update 
```

* poreCov with a basecalled fastq directory 

```bash
nextflow run replikation/poreCov --fastq_pass 'fastq_pass/' -r 0.8.0 \
    --cores 32  -profile local,docker
```

* poreCov with basecalling and renaming of barcodes based on `sample_names.csv`

```bash
# rename barcodes automatically by providing an input file, also using another primer scheme
nextflow run replikation/poreCov --fast5 fast5_dir/ --samples sample_names.csv \
   --primerV V1200 --output results -profile local,docker
```

## 2.3 Extended Usage
* see also `nextflow run replikation/poreCov --help -r 0.8.0`
### Version control
* poreCov supports version control via `-r` this way, you can run everything reproducible (e.g. `-r 0.8.0`)
* poreCov releases are listed [here](https://github.com/replikation/poreCov/releases)
* add `-r <version>` to a poreCoV run to activate this
* run `nextflow pull replikation/poreCov` to install updates
   * if you have issues during update try `rm -rf ~/.nextflow` and then `nextflow pull replikation/poreCov`
   * this removes old files and downloads everything new

### Important input flags (choose on)
* these are the flags to get "data" into the workflow
   * `--fast5 fast5_dir/`  for fast5 directory input
   * `--fastq_pass fastq_dir/`  directory with basecalled data (contains "barcode01" etc. directories)
   * `--fastq "sample*.fastq.gz"` alternative fastq input (one sample per file)
   * `--fasta "*genomes.fasta"`  SARS-CoV-2 genomes as fasta 

### Sample sheet
* barcodes can be automatically renamed via `--samples sample_names.csv`
* example comma separated file (don't replace the header)
  * `_id` = sample name
  * `Status` = barcode number which should be renamed

```csv
_id,Status
Sample_2021,barcode01
2ndSample,BC02
```
  
### Pangolin Lineage definitions
  * lineage determinations are quickly changing in response to the pandemic
  * to avoid using out of date lineage schemes, a `--update` flag can be added to each poreCov run to get the most recent version-controlled pangolin container
  * we are currently building two times every week version-controlled pangolin container automatically, [see here](https://hub.docker.com/r/nanozoo/pangolin/tags?page=1&ordering=last_updated)
    * it is also possible to use instead of `--update` this flag: `--pangolindocker  'nanozoo/pangolin:3.1.1--2021-06-14'`
    * this way you can use other container or version, but beware some containers might not be compatible with poreCov
  

# 3. Quality Metrics (default)

* Regions with coverage of 20 or less are masked ("N")
* Genome quality is compared to NC_045512.2
    * Genome quality assessment is based on [RKIBioinformaticsPipelines/president](https://gitlab.com/RKIBioinformaticsPipelines/president)
        * also prepares csv and fasta for upload via DESH portal
* Pangolin lineages are determined
* nextstrain clades are determined including mutation infos
* reads are classified to human and SARS-CoV-2 to check for possible contamination and sample prep issues


# 4. Workflow

* poreCov was coded with "easy to use" in mind, while staying flexible
* therefore we provide a few input types which adjusts the workflow automatically (see image below)
  * fast5 raw data, fastq files (one sample per file), fastq_pass (the basecalling output) or fasta (supports multifastas)
* primer schemes for ARTIC can be V1, V2, V3(default) or V1200 (the 1200bp amplicon ones)

<p align="left">
  <img src="data/figures/workflow.png" width="700" title="Workflow">
</p>


# 5. Literature / References to cite
For citing etc. check out these programs used for poreCov:
* [artic protocol](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)
* [kraken2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)
* [krona](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-385)
* [medaka](https://github.com/nanoporetech/medaka)
* [minimap2](https://github.com/lh3/minimap2)
* [nextclade](https://clades.nextstrain.org/)
* [Nextflow](https://www.nextflow.io/index.html)
* [pangolin](https://github.com/hCoV-2019/pangolin)
* [president](https://gitlab.com/RKIBioinformaticsPipelines/president)
* [CoVarPlot](https://github.com/Psy-Fer/CoVarPlot)

# 6. Troubleshooting
* Collection of some helpful infos

## Singularity
* Singularity needs additional option flags to run like `--userns` [Solution on how to pass Singularity commands to poreCov](https://github.com/replikation/poreCov/issues/101#issuecomment-825807042)
  
# 7. Time to results
**Table 1**: Execution speed of poreCov on different Ubuntu 20 Systems using a single sample file with 167,929 reads. Command used: `nextflow run replikation/poreCov -profile test_fastq,local,docker`.
Hardware|First time with download (DB+container)¹ |Default settings | 
-|-|-
2 CPUs 4 GB RAM| 1h 2min | 32 min 30s ² 
2 CPUs 8 GB RAM| 46 min | 21m 20s 
4 CPUs 16 GB RAM| 40 min | 12m 48s 
8 CPUs 32  GB RAM| 35 min | 11m 39s 
16 CPUs 64  GB RAM| 30 min | 9m 39s

¹ *time depends mostly on available internet speed*  
² *was not able to execute read classification due to limited hardware, but generated and classified SARS-CoV-2 genomes*  


**Table 2**: Execution speed of poreCov on different Ubuntu 20 Systems using 24 fastq samples. Command used: `nextflow run replikation/poreCov -profile test_fastq,local,docker`. Time meassured by the start of the workflow.
Hardware|Default settings 
-|-|
2 CPUs 4 GB RAM| 13h 33m ¹
2 CPUs 8 GB RAM| 7h 56m  ¹
4 CPUs 16 GB RAM| 4h 10 min
8 CPUs 32  GB RAM| 2h 15 min
16 CPUs 64  GB RAM| 1h 25 min
  
¹ *was not able to execute read classification due to limited hardware, but generated and classified SARS-CoV-2 genomes*  
