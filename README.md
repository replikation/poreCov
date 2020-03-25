![logo](data/logo/mobile_logo.png)
**nCov Workflows**   |   Institut für Infektionsmedizin und Krankenhaushygiene
===

![](https://img.shields.io/badge/nextflow-20.01.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)
![](https://github.com/replikation/nCov/workflows/Syntax_check/badge.svg)


[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 

* by Christian Brandt
* **this tool is currently under heavy development, so expect some bugs but feel free to report issues**

## What is this Repo?

* an attempt to streamline the usage of various nCov analysis workflows
* the main focus is stability and easy to use for the User

## Installation

* nCov runs with the workflow manager `nextflow` using `docker`
* this means all the other programs are automatically pulled via docker
* Only `docker` and `nextflow` needs to be installed

### Easy Installation
* if you dont have experience with bioinformatic tools use this
* just copy the commands into your terminal to set everything up

```bash
sudo apt-get update
sudo apt install -y default-jre
curl -s https://get.nextflow.io | bash 
sudo mv nextflow /bin/
sudo apt-get install -y docker-ce docker-ce-cli containerd.io
sudo usermod -a -G docker $USER
```

* restart your computer

* try out the installation by entering the following

```bash
nextflow run replikation/nCov --help
```

### Normal Installation

* this is the default choice

**Dependencies**

>   * docker (add docker to your Usergroup, so no sudo is needed)
>   * nextflow + java runtime 
>   * git (should be already installed)
>   * wget (should be already installed)
>   * tar (should be already installed)

* Docker installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce)
* Nextflow installation [here](https://www.nextflow.io/)
* move or add the nextflow executable to a bin path
* add docker to your User group via `sudo usermod -a -G docker $USER`
