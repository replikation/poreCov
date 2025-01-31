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

[![](https://img.shields.io/badge/publication-frontiers-brightgreen.svg)](https://doi.org/10.3389/fgene.2021.711437)

[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 


**Citation:**  
> **poreCov - an easy to use, fast, and robust workflow for SARS-CoV-2 genome reconstruction via nanopore sequencing**  
> Christian Brandt, Sebastian Krautwurst, Riccardo Spott, Mara Lohde, Mateusz Jundzill, Mike Marquet, Martin Hölzer  
> https://www.frontiersin.org/articles/10.3389/fgene.2021.711437/full  

## What is this Repo?

* poreCov is a SARS-CoV-2 analysis workflow for nanopore data (via the [ARTIC protocol](https://artic.network/ncov-2019)) or SARS-CoV-2 genomes (fasta)
* the workflow is pre-configured to simplify [data analysis](https://htmlpreview.github.io/?https://github.com/replikation/poreCov/blob/master/data/figures/index.html): 
<p align="left">
    <a href="https://htmlpreview.github.io/?https://github.com/replikation/poreCov/blob/master/data/figures/index.html">
        <img src="data/figures/report_summary.png" width="500" title="Report file">
</p>

Table of Contents
=================
<!--ts-->
- [**poreCov | SARS-CoV-2 Workflow for nanopore sequencing data**](#porecov--sars-cov-2-workflow-for-nanopore-sequencing-data)
  - [What is this Repo?](#what-is-this-repo)
- [Table of Contents](#table-of-contents)
- [1. Quick Setup (Ubuntu)](#1-quick-setup-ubuntu)
  - [1.1 Nextflow (the workflow manager)](#11-nextflow-the-workflow-manager)
  - [1.2 Container (choose one - they manage all the tools)](#12-container-choose-one---they-manage-all-the-tools)
    - [Docker](#docker)
    - [Singularity](#singularity)
    - [Conda (not recommended)](#conda-not-recommended)
  - [1.3 Basecalling (optional)](#13-basecalling-optional)
- [2. Run poreCov](#2-run-porecov)
  - [2.1 Test run](#21-test-run)
  - [2.2 Quick run examples](#22-quick-run-examples)
  - [2.3 Extended Usage](#23-extended-usage)
    - [Version control](#version-control)
    - [Important input flags (choose one)](#important-input-flags-choose-one)
    - [Custom primer bed files](#custom-primer-bed-files)
    - [Sample input](#sample-input)
      - [Sample sheet](#sample-sheet)
      - [List input](#list-input)
    - [Pangolin Lineage definitions](#pangolin-lineage-definitions)
- [3. Quality Metrics (default)](#3-quality-metrics-default)
- [4. Workflow](#4-workflow)
- [5. Literature / References to cite](#5-literature--references-to-cite)
- [6. Troubleshooting](#6-troubleshooting)
  - [Singularity](#singularity-1)
- [7. Time to results](#7-time-to-results)
- [8. Credits](#8-credits)
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
* Singularity installation [here](https://apptainer.org/docs/)
* if you can't use Docker

Note, that with Singularity the following environment variables are automatically passed to the container to ensure execution on HPCs: `HTTPS_PROXY`, `HTTP_PROXY`, `http_proxy`, `https_proxy`, `FTP_PROXY` and `ftp_proxy`.
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
nextflow run replikation/poreCov -profile test_fastq,local,docker -r 1.1.0 --update

# or for Singularity or conda installation
nextflow run replikation/poreCov -profile test_fastq,local,singularity -r 1.1.0 --update
```

## 2.2 Quick run examples

* poreCov with basecalling and Docker
    * `--update` tryies to force the most recent pangolin lineage and nextclade release version (optional)
    * `-r 1.1.0` specifies the workflow release from [here](https://github.com/replikation/poreCov/releases)
    * `--primerV` specifies the primer sets that were used, see `--help` to see what is supported
        * alternatively provide a primer bed file on your own
```bash
nextflow run replikation/poreCov --fast5 fast5/ -r 1.1.0 \
    --cores 6 -profile local,docker --update --primerV V4
```

* poreCov with a basecalled fastq directory and custom primer bed file

```bash
nextflow run replikation/poreCov --fastq_pass 'fastq_pass/' -r 1.1.0 \
    --cores 32  -profile local,docker --update --primerV primers.bed
```

* poreCov with basecalling and renaming of barcodes based on `sample_names.csv`

```bash
# rename barcodes automatically by providing an input file, also using another primer scheme
nextflow run replikation/poreCov --fast5 fast5_dir/ --samples sample_names.csv \
   --primerV V1200 --output results -profile local,docker --update
```

## 2.3 Extended Usage
* see also `nextflow run replikation/poreCov --help -r 1.1.0`
### Version control
* poreCov supports version control via `-r` this way, you can run everything reproducible (e.g. `-r 1.1.0`)
    * moreover only releases are extensively tested and validated
* poreCov releases are listed [here](https://github.com/replikation/poreCov/releases)
* add `-r <version>` to a poreCoV run to activate this
* run `nextflow pull replikation/poreCov` to install updates
   * if you have issues during update try `rm -rf ~/.nextflow` and then `nextflow pull replikation/poreCov`
   * this removes old files and downloads everything new

### Important input flags (choose one)
* these are the flags to get "data" into the workflow
   * `--fast5 fast5_dir/`  for fast5 directory input
   * `--fastq_pass fastq_dir/`  directory with basecalled data (contains "barcode01" etc. directories)
   * `--fastq "sample*.fastq.gz"` alternative fastq input (one sample per file)
   * `--fasta "*genomes.fasta"`  SARS-CoV-2 genomes as fasta (.gz allowed)

### Custom primer bed files
* poreCov supports the input of `primer.bed` files via `--primerV` instead of selecting a preexisting primer version like `--primerV V4`
   * for an example see [2.2 Quick run examples](#22-quick-run-examples)
   * feature available for poreCov version `1.1.0` or greater
* the main issue with primer bed files is that they need to have the correct columns and text to be recognized via artic
* the following rules apply to the bed file (see also example)
    * each column is separated via one `tab` or `\t`
    * column 1 is the fasta reference, and it should be MN908947.3 (poreCov replaces that automatically)
    * column 2 is the primer start
    * column 3 is the primer end
    * column 4 is the primer name, and it *has to end* with `_RIGHT` or `_LEFT`
    * column 5 is the pool and it should be named `nCoV-2019_1` or `nCoV-2019_2`
    * column 6 defines the strand orientation with either `-` or `+`

```csv
MN908947.3	30	54	nCoV-2019_1_LEFT	nCoV-2019_1	+
MN908947.3	1183	1205	nCoV-2019_1_RIGHT	nCoV-2019_1	-
MN908947.3	1100	1128	nCoV-2019_2_LEFT	nCoV-2019_2	+
MN908947.3	2244	2266	nCoV-2019_2_RIGHT	nCoV-2019_2	-
MN908947.3	2153	2179	nCoV-2019_3_LEFT	nCoV-2019_1	+
MN908947.3	3235	3257	nCoV-2019_3_RIGHT	nCoV-2019_1	-
MN908947.3	3144	3166	nCoV-2019_4_LEFT	nCoV-2019_2	+
MN908947.3	4240	4262	nCoV-2019_4_RIGHT	nCoV-2019_2	-
```

### Sample input

> [!NOTE]  
> If using --fastq without either --sample or --list, samples whose concatenated and size-selected FastQ files are smaller than 1500 kB will be excluded from further analysis.

#### Sample sheet
* barcodes can be automatically renamed via `--samples sample_names.csv`
* required columns:
  * `_id` = sample name
  * `Status` = barcode number which should be renamed
* optional column:
  * `Description` = description column to be included in the output report and tables

Example comma separated file (don't replace the header):
```csv
_id,Status,Description
Sample_2021,barcode01,good
2ndSample,BC02,bad
```

#### List input
* Using `--list` You can provide a csv as input to `--fastq` to select for specific fastq-files
  * e.g.: `--fastq input.csv --list`
  * the csv needs to contain two columns:
    * column 1 = sample name
    * column 2 = path to fastq-location
  * no header should be used
* files get automatically renamed to the sample names provided in column 1

 Example:
```csv
sample1,path/to/first/sample.fastq.gz
2ndSample,path/to/second/sample.fastq.gz
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
* primer schemes for ARTIC can be V1, V2, V3(default), V4, V4.1 (the 400bp amplicon ones), V1200 or V5.2.0_1200 (the 1200bp amplicon ones)

<p align="left">
  <img src="data/figures/workflow.png" width="700" title="Workflow">
</p>


# 5. Literature / References to cite
If you are using poreCov please also check the used software to cite in your work:
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
* [LCS](https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btac047/42357424/btac047.pdf)

# 6. Troubleshooting
* Collection of some helpful infos

## Singularity
* Singularity needs additional option flags to run like `--userns` [Solution on how to pass Singularity commands to poreCov](https://github.com/replikation/poreCov/issues/101#issuecomment-825807042)
  
# 7. Time to results
**Table 1**: Execution speed of poreCov on different Ubuntu 20 Systems using a single sample file with 167,929 reads. Command used: `nextflow run replikation/poreCov -r 0.9.4 -profile test_fastq,local,docker`.
  
Hardware|First time with download (DB+container)¹ |Default settings | 
-|-|-
2 CPUs 4 GB RAM| 1h 2min | 32 min 30s ² 
2 CPUs 8 GB RAM| 46 min | 21m 20s 
4 CPUs 16 GB RAM| 40 min | 12m 48s 
8 CPUs 32  GB RAM| 35 min | 11m 39s 
16 CPUs 64  GB RAM| 30 min | 9m 39s

¹ *time depends mostly on available internet speed*  
² *was not able to execute read classification due to limited hardware, but generated and classified SARS-CoV-2 genomes*  


**Table 2**: Execution speed of poreCov on different Ubuntu 20 Systems using 24 fastq samples. Command used: `nextflow run replikation/poreCov -r 0.9.4 --fastq "*.fastq.gz" --primerV V1200 --samples samplenames.csv -profile local,docker`. Time meassured by the start of the workflow.
  
Hardware|Default settings 
-|-|
2 CPUs 4 GB RAM| 13h 33m ¹
2 CPUs 8 GB RAM| 7h 56m  ¹
4 CPUs 16 GB RAM| 4h 10 min
8 CPUs 32  GB RAM| 2h 15 min
16 CPUs 64  GB RAM| 1h 25 min
  
¹ *was not able to execute read classification due to limited hardware, but generated and classified SARS-CoV-2 genomes*  
  
# 8. Credits
The key steps of poreCov are carried out using the [ARTIC Network field bioinformatics pipeline](https://github.com/artic-network/fieldbioinformatics). Kudos to all amazing developers for your incredible efforts during this pandemic! Many thanks to all others who have helped out and contributed to poreCov as well. 

The script [convert_VCF_info_fields.py](bin/convert_VCF_info_fields.py) was originally developed by the Intergalactic Utilities Commission in the [Tool Shed repository](https://github.com/galaxyproject/tools-iuc) under a MIT license.