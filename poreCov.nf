#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* Nextflow -- SARS-CoV-2 Analysis Pipeline
* Author: christian.jena@gmail.com
*/

/************************** 
* HELP messages & checks
**************************/

/* 
Nextflow version check  
Format is this: XX.YY.ZZ  (e.g. 20.07.1)
change below
*/

XX = "20"
YY = "07"
ZZ = "1"

if ( nextflow.version.toString().tokenize('.')[0].toInteger() < XX.toInteger() ) {
println "\033[0;33mporeCov requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}
else if ( nextflow.version.toString().tokenize('.')[1].toInteger() == XX.toInteger() && nextflow.version.toString().tokenize('.')[1].toInteger() < YY.toInteger() ) {
println "\033[0;33mporeCov requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}

// Log infos based on user inputs
if ( params.help ) { exit 0, helpMSG() }
    defaultMSG()
if ( params.primerV.matches('V1200') ) { v1200_MSG() }
if ( params.dir || workflow.profile.contains('test_fast5') ) { basecalling() }

// profile helps
    if ( workflow.profile == 'standard' ) { exit 1, "NO EXECUTION PROFILE SELECTED, use e.g. [-profile local,docker]" }
    if (params.profile) { exit 1, "--profile is WRONG use -profile" }
    if (
        workflow.profile.contains('singularity') ||
        workflow.profile.contains('nanozoo') ||
        workflow.profile.contains('ukj_cloud') ||
        workflow.profile.contains('docker')
        ) { "engine selected" }
    else { println "No engine selected:  -profile EXECUTER,ENGINE" 
           println "using native installations" }
    if (
        workflow.profile.contains('nanozoo') ||
        workflow.profile.contains('ukj_cloud') ||
        workflow.profile.contains('local')
        ) { "executer selected" }
    else { exit 1, "No executer selected:  -profile EXECUTER,ENGINE" }

    if (workflow.profile.contains('local')) {
        println "\033[2m Using $params.cores/$params.max_cores CPU threads [--max_cores]\u001B[0m"
        println " "
    }
    if ( workflow.profile.contains('singularity') ) {
        println ""
        println "\033[0;33mWARNING: Singularity image building sometimes fails!"
        println "Multiple resumes (-resume) and --max_cores 1 --cores 1 for local execution might help.\033[0m\n"
    }

// params help
if (!workflow.profile.contains('test_fastq') && !workflow.profile.contains('test_fast5') && !workflow.profile.contains('test_fasta')) {
    if (!params.fasta &&  !params.dir &&  !params.fastq ) {
        exit 1, "input missing, use [--fasta] [--fastq] or [--dir]"}
    if ((params.fasta && ( params.fastq || params.dir )) || ( params.fastq && params.dir )) {
        exit 1, "To much inputs: please us either: [--fasta], [--fastq] or [--dir]"} 
    if (!params.metadata) { println "\033[0;33mNo [--metadata] file specified, skipping tree build\u001B[0m" }
}
/************************** 
* INPUTs
**************************/

// fasta input 
    if (params.fasta && !workflow.profile.contains('test_fasta')) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
    }

// consensus qc reference input - auto using git default if not specified
    if (params.reference_for_qc) { 
        reference_for_qc_input_ch = Channel
        .fromPath( params.reference_for_qc, checkIfExists: true)
    }
    else if (!params.reference_for_qc) {
        reference_for_qc_input_ch = Channel
        .fromPath(workflow.projectDir + "/data/reference_nCov19/NC_045512.2.fasta")
    }

// references input 
    if (params.references) { reference_input_ch = Channel
        .fromPath( params.references, checkIfExists: true)
    }

// metadata input 
    if (params.metadata) { metadata_input_ch = Channel
        .fromPath( params.metadata, checkIfExists: true)
    }

// fastq input or via csv file
    if (params.fastq && params.list && !workflow.profile.contains('test_fastq')) { 
        fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
    }
    else if (params.fastq && !workflow.profile.contains('test_fastq')) { 
        fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
    }

// dir input
    if (params.dir && !workflow.profile.contains('test_fast5')) { dir_input_ch = Channel
        .fromPath( params.dir, checkIfExists: true, type: 'dir')
        .map { file -> tuple(file.name, file) }
    }

/************************** 
* MODULES
**************************/

include { get_nanopore_fastq } from './modules/get_fastq_test_data.nf'
include { get_fasta } from './modules/get_fasta_test_data.nf'
include { get_fast5 } from './modules/get_fast5_test_data.nf'

/************************** 
* Workflows
**************************/

include { artic_ncov_wf } from './workflows/artic_nanopore_nCov19.nf'
include { basecalling_wf } from './workflows/basecalling.nf'
include { build_database_wf } from './workflows/databases.nf'
include { create_tree_wf } from './workflows/create_tree.nf'
include { determine_lineage_wf } from './workflows/determine_lineage.nf'
include { genome_quality_wf } from './workflows/genome_quality.nf'
include { read_qc_wf } from './workflows/read_qc.nf'
include { toytree_wf } from './workflows/toytree.nf'

/************************** 
* MAIN WORKFLOW
**************************/

workflow {
    // 0. Test profile data
        if ( workflow.profile.contains('test_fastq')) { fastq_input_ch =  get_nanopore_fastq().map {it -> ['SARSCoV2', it] } }
        if ( workflow.profile.contains('test_fasta')) { fasta_input_ch =  get_fasta().map {it -> ['SARSCoV2', it] } }
        if ( workflow.profile.contains('test_fast5')) { dir_input_ch =  get_fast5().map {it -> ['SARSCoV2', it] } }

    // 1. Reconstruct genomes
        if (params.dir || workflow.profile.contains('test_fast5')) { 
            fasta_input_ch = artic_ncov_wf(basecalling_wf(dir_input_ch))
        }
        if (params.fastq || workflow.profile.contains('test_fastq')) { 
            read_qc_wf(fastq_input_ch)
            fasta_input_ch = artic_ncov_wf(fastq_input_ch)
        }

    // 2. Genome quality and lineages
        determine_lineage_wf(fasta_input_ch)
        genome_quality_wf(fasta_input_ch, reference_for_qc_input_ch)


    // 3. (optional) analyse genomes to references and build tree
        if (params.references && params.metadata && (params.fastq || params.fasta || params.dir)) {
        // build tree 
            create_tree_wf (fasta_input_ch, reference_input_ch, metadata_input_ch) 
                newick = create_tree_wf.out
        }

        else if (params.metadata && (params.fastq || params.fasta || params.dir)) {
        // build database
            build_database_wf()
        // merge build_database_wf metadata with user metadata file
            meta_merge_ch = build_database_wf.out[1].splitCsv(header: true, sep: '\t')
                .mix(metadata_input_ch.splitCsv(header: true, sep: '\t'))
                .collectFile(seed: 'strain\tcountry\tdate\n') { 
                    row -> [ "metadata.tsv", row.strain + '\t' + row.country + '\t' + row.date + '\n' ]  }
        // build tree
            create_tree_wf (fasta_input_ch, build_database_wf.out[0], meta_merge_ch)
                newick = create_tree_wf.out
            toytree_wf(newick)
        }
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
    
    ${c_green}poreCov${c_reset} | A Nextflow SARS-CoV-2 (nCov19) workflow for nanopore data
    
    ${c_yellow}Usage examples:${c_reset}
    nextflow run replikation/poreCov --fastq 'sample_01.fasta.gz' --cores 14 -profile local,singularity

    ${c_yellow}Inputs (choose one):${c_reset}
    --dir           one fast5 dir of a nanopore run containing multiple samples (barcoded);
                    to run single sample (no barcodes) add [--single]
                    ${c_dim}(not implemented yet) autorename barcodes via [--barcodeIDs rename.csv] 
                        Per line: 01;samplename${c_reset}
                    ${c_dim}[basecalling - demultiplexing - nCov genome reconstruction]${c_reset}

    --fastq         one fastq or fastq.gz file per sample or
                    multiple file-samples: --fastq 'sample_*.fasta.gz'
                    ${c_dim}[nCov genome reconstruction]${c_reset}

    --fasta         direct input of genomes, one file per genome
                    ${c_dim}[Lineage determination, Quality control]${c_reset}

    ${c_yellow}Parameters - Basecalling${c_reset}
    --localguppy    use a native guppy installation instead of a gpu-guppy-docker 
                    native guppy installation is used by default for singularity or conda
    --one_end       removes the recommended "--require_barcodes_both_ends" from guppy demultiplexing
                    try this if to many barcodes are unclassified (check the pycoQC report)
    --guppy_cpu     use cpus instead of gpus for basecalling

    ${c_yellow}Parameters - nCov genome reconstruction${c_reset}
    --primerV       artic-ncov2019 primer_schemes [default: ${params.primerV}]
                        Supported: V1, V2, V3, V1200
    --minLength     min length filter raw reads [default: ${params.minLength}]
    --maxLength     max length filter raw reads [default: ${params.maxLength}]

    ${c_yellow}Parameters - Genome quality control${c_reset}
    --reference_for_qc      reference FASTA for consensus qc (optional, wuhan is provided by default)
    --threshold             global pairwise sequence identity threshold [default: ${params.threshold}] 

    ${c_yellow}Parameters - Tree construction:${c_reset}
    Input is either: --fasta --fastq --dir

    --references    multifasta file to compare against your input

    --metadata      tsv file with 3 rows and header: strain country date   
                    date in YYYY-MM-DD   strain is fasta header without >
    
    Optional:
    --highlight     names containing this string are colored in the tree in red 
                    [default: ${params.highlight}]
    --maskBegin     masks beginning of alignment [default: ${params.maskBegin}]
    --maskEnd       masks end of alignment [default: ${params.maskEnd}]
    --rm_N_genome   removes genomes from tree with x amount of N's or more [default: ${params.rm_N_genome}]

    ${c_yellow}Options:${c_reset}
    --cores         max cores for local use [default: $params.cores]
    --max_cores     max cores used on the machine for local use [default: $params.max_cores]
    --memory        available memory [default: $params.memory]
    --output        name of the result folder [default: $params.output]
    --cachedir      defines the path where singularity images are cached
                    [default: $params.cachedir] 

    ${c_yellow}Execution/Engine profiles:${c_reset}
    poreCov supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.:
     -profile ${c_green}local${c_reset},${c_blue}docker${c_reset}
     -profile ${c_yellow}test_fastq${c_reset},${c_green}slurm${c_reset},${c_blue}singularity${c_reset}

      ${c_green}Executer${c_reset} (choose one):
       local
       slurm
      ${c_blue}Engines${c_reset} (choose one):
       docker
       singularity
      ${c_yellow}Input test data${c_reset} (choose one):
       test_fasta
       test_fastq
       test_fast5

    Alternatively provide your own configuration via -c ownconfig.config 
    """.stripIndent()
}

def defaultMSG(){
    log.info """
    SARS-CoV-2 - Workflow


    \u001B[32mProfile:             $workflow.profile\033[0m
    \033[2mCurrent User:        $workflow.userName
    Nextflow-version:    $nextflow.version
    \u001B[0m
    Pathing:
    \033[2mWorkdir location [-work-Dir]:
        $workflow.workDir
    Output dir [--output]: 
        $params.output
    Singularity cache dir [--cachedir]: 
        $params.cachedir
    \u001B[1;30m______________________________________\033[0m
    Parameters:
    \033[2mPrimerscheme:        $params.primerV [--primerV]
    CPUs to use:         $params.cores [--cores]
    Memory in GB:        $params.memory [--memory]\u001B[0m
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}

def v1200_MSG() {
    log.info """
    1200 bp options are used as primer scheme (V1200)
      --minLength set to 250bp
      --maxLength set to 1500bp
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}

def basecalling() {
    log.info """
    Basecalling options:
    \033[2mUsing local guppy?      $params.localguppy [--localguppy]  
    One end demultiplexing? $params.one_end [--one_end]
    CPUs for basecalling?   $params.guppy_cpu [--guppy_cpu]\u001B[0m
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}
