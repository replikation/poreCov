#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* Nextflow -- nCov Analysis Pipeline
* Author: christian.jena@gmail.com
*/

/************************** 
* HELP messages & USER INPUT checks
**************************/
if( !nextflow.version.matches('20.+') ) {
    println "This workflow requires Nextflow version 20.X or greater -- You are running version $nextflow.version"
    exit 1
}

if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mnCov - Workflows\033[0m"
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
println "\033[2mCPUs to use: $params.cores"
println "\033[2mMemory in GB: $params.memory"
if (params.dir) { println "\033[2mBarcodes: $params.barcodes" }
if (params.artic_ncov19) { println "\033[2mPrimerscheme: $params.primerV"  }
println "Output dir: $params.output\u001B[0m"
println " "

if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (!params.fasta &&  !params.dir &&  !params.fastq ) {
    exit 1, "input missing, use [--fasta] [--fastq] or [--dir]"}
if (params.fasta && ( params.fastq || params.dir )) {
    exit 1, "To much inputs: please us either: [--fasta], [--fastq] or [--dir]"} 
if (params.augur && (!params.metadata || !params.references)) {
    exit 1, "Please provide for augur: [--references] and [--metadata]"} 

// fasta input 
    if (params.fasta) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    }

// references input 
    if (params.references) { reference_input_ch = Channel
        .fromPath( params.references, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
    }

// metadata input 
    if (params.metadata) { metadata_input_ch = Channel
        .fromPath( params.metadata, checkIfExists: true)
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
include augur_align from './modules/augur'
include augur_tree from './modules/augur'
include augur_tree_refine from './modules/augur'
include fasttree from "./modules/fasttree"
include filter_fastq_by_length from './modules/filter_fastq_by_length'
include guppy_gpu from './modules/guppy'
include mafft from './modules/mafft'
include mask_alignment from './modules/mask_alignment'
include snp_sites from './modules/snp_sites'
include split_reference from './modules/split_reference'
include toytree from './modules/toytree'

/************************** 
* SUB WORKFLOWS
**************************/



workflow dir_handler_wf {
    take: 
        dir_input  
    main:
        
        guppy_gpu(dir_input)
        
        if (params.barcodes) {
            fastq_channel = guppy_gpu.out
                            .map { it -> it[1] }
                            .flatten()
                            .map { it -> [ it.simpleName, it ] }
            }
        else { fastq_channel = guppy_gpu.out }
    
    emit:
        fastq_channel
} 


workflow artic_nCov19_wf {
    take:   
        fastq
    main:   
        artic(filter_fastq_by_length(fastq))
    emit:   
        artic.out
}

workflow create_tree_nextstrain_wf {
    take: 
        fasta       // the nCov fasta (own samples or reconstructed here)
        references  // multiple references to compare against
        metadata    // tsv file of meta data  strain country date
    main:

        align_reference = Channel.fromPath( workflow.projectDir + "/data/reference_nCov19/MN908947.gb", checkIfExists: true)

        augur_tree(
            mask_alignment(
                augur_align(fasta, references, align_reference)))

        augur_tree_refine(augur_tree.out, metadata)

    emit:
        augur_tree_refine.out
}

workflow create_tree_mafft_wf {
    take: 
        fasta
        references
    main:
        fasttree(
                mafft (fasta.map{it -> it [1] }.collect(), references))
            
    emit:
        fasttree.out
}

// TODO: get fastaname and carry it as env / val to the toytree highlight
// this way i can highlight all samples in there

// also you could parse a highlight syntax so ppl can do --highlight "UKJ" therefore annotate multiple key words
workflow toytree_wf {
    take: 
        trees  
    main:
        toytree(trees.map{ it -> [it[0], it[2]] })
    emit:
        toytree.out
} 

/************************** 
* MAIN WORKFLOW
**************************/

workflow {


// get genome workflows
    if (params.artic_ncov19 && params.dir) { artic_nCov19_wf(dir_handler_wf(dir_input_ch)); fasta_input_ch = artic_nCov19_wf.out }
    if (params.artic_ncov19 && params.fastq) { artic_nCov19_wf(fastq_input_ch); fasta_input_ch = artic_nCov19_wf.out}

// analyse genome to references
    if (params.references && (params.fastq || params.fasta || params.dir)) { 
        if (params.mafft) { create_tree_mafft_wf (fasta_input_ch, reference_input_ch); newick = create_tree_mafft_wf.out } 
        if (params.augur)  { create_tree_nextstrain_wf (fasta_input_ch, reference_input_ch, metadata_input_ch); newick = create_tree_nextstrain_wf.out }

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
    
    Nextflow nCov workflows, by Christian Brandt
    
    ${c_yellow}Usage examples:${c_reset}
    nextflow run replikation/nCov --artic_ncov19 --fastq 'sample_01.fasta.gz'

    nCov.nf --artic_ncov19 --fastq 'sample_01.fasta.gz' \\ 
                    --augur --references references.fasta \\
                    --metadata metadata.tsv \\
                    --cores 8 -profile local,docker

    ${c_yellow}Reconstruct genome workflows:${c_reset}
    ${c_blue} --artic_ncov19 ${c_reset}
    ${c_dim}Inputs:${c_reset}
    ${c_dim}--fastq${c_reset}   one fastq file per sample as input
    ${c_dim}--dir     one fast5 dir of a nanopore run 
                                add [--barcodes] to the command if they are barcoded${c_reset}
                      
    ${c_dim}Parameters:${c_reset} 
    ${c_dim} --primerV       artic-ncov2019 primer_schemes [default: ${params.primerV}]${c_reset}
    ${c_dim} --minLength     min length filter raw reads [default: ${params.minLength}]${c_reset}
    ${c_dim} --maxLength     max length filter raw reads [default: ${params.maxLength}]${c_reset}
 
    ${c_yellow}Pyholgenetic tree workflows:${c_reset}
    ${c_blue} --augur${c_reset}         Input own genomes via ${c_green}[--fasta]${c_reset}
                     or use ${c_green}[--fastq],[--dir]${c_reset} with  ${c_blue}--<reconstruct workflow>${c_reset}
    ${c_dim}Mandatory: 
    ${c_dim} --references    fasta file(s) to compare against your data sample 
    ${c_dim} --metadata      tsv file with 3 headers: strain country date   (date in YYYY-MM-DD)
    ${c_dim}Parameters: 
    ${c_dim} --maskBegin     masks beginning of alignment [default: ${params.maskBegin}]
    ${c_dim} --maskEnd       masks end of alignment [default: ${params.maskEnd}]

    ${c_blue} --mafft${c_reset}         Input own genomes via ${c_green}[--fasta]${c_reset} or 
                     use ${c_green}[--fastq],[--dir]${c_reset} with  ${c_blue}--<reconstruct workflow>${c_reset}

    ${c_reset}Options:
    --cores          max cores for local use [default: $params.cores]
    --memory         available memory [default: $params.memory]
    --output         name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html       cpu / ram usage (may cause errors)
    -with-dag chart.html        generates a flowchart for the process tree
    -with-timeline time.html    timeline (may cause errors)

    Profile:
    -profile                 local,docker -> merge profiles e.g. -profile local,docker ${c_reset}
    """.stripIndent()
}