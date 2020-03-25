#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* Nextflow -- nCov Analysis Pipeline
* Author: christian.jena@gmail.com
*/

/************************** 
* HELP messages & USER INPUT checks
**************************/
if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (!params.fasta &&  !params.dir &&  !params.fastq ) {
    exit 1, "input missing, use [--fasta] [--fastq] or [--dir]"}
if (params.fasta && params.fastq) {
    exit 1, "please us either: [--fasta] or [--fastq]"}   

// fasta input 
    if (params.fasta) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    }

// fastq input
    if (params.fastq) { fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    }

// dir input
    if (params.dir) { dir_input_ch = Channel
        .fromPath( params.dir, checkIfExists: true, type: 'dir')
        .map { file -> tuple(file.name, file) }
    }


/************************** 
* DATABASES
**************************/






/************************** 
* MODULES
**************************/
    include artic from './modules/artic' 
    include filter_fastq_by_length from './modules/filter_fastq_by_length'
    include cat_fastq from './modules/cat_fastq'

/************************** 
* SUB WORKFLOWS
**************************/

workflow artic_nCov19_wf {
    take:   fastq
    main:   artic(filter_fastq_by_length(fastq))
    emit:   medaka.out
}


/************************** 
* MAIN WORKFLOW
**************************/

workflow {
    if (params.artic_ncov19 && params.dir) { artic_nCov19_wf(cat_fastq(dir_input_ch)) fasta_input_ch = artic_nCov19_wf.out }
    if (params.artic_ncov19 && params.fastq) { artic_nCov19_wf(fastq_input_ch); fasta_input_ch = artic_nCov19_wf.out}

    if (params.references && (params.fastq || params.fasta)) { build_tree_wf }
}

/*************  
* --help
*************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    Nextflow nCov workflows for easy use, by Christian Brandt
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run replikation/nCov ${c_blue}--artic_ncov19${c_reset} ${c_green}--fastq 'sample_01.fasta.gz'${c_reset} -profile local,docker

    ${c_yellow}Workflow options:${c_reset}
    ${c_blue} --artic_ncov19 ${c_reset}
    Input options:
    ${c_green}--fastq Sample_01.fastq${c_reset}     (one sample per fastq file)
    ${c_green}--fastq 'Sample_*.fastq'${c_reset}    (multiple samples at once, one sample per fastq)
    ${c_green}--dir fastq_files/${c_reset}          (on dir containing multiple fastq files for one sample)
    Setting:
    ${c_green}--primerV ${params.primerV}${c_reset} (artic-ncov2019 primer_schemes version used, available: V1, V2, V3)


    ${c_reset}Options:
    --cores             max cores for local use [default: $params.cores]
    --output            name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 local,docker -> merge profiles e.g. -profile local,docker ${c_reset}
    """.stripIndent()
}