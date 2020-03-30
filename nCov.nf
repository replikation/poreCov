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

// references input 
    if (params.references) { reference_input_ch = Channel
        .fromPath( params.references, checkIfExists: true)
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
include cat_fastq from './modules/cat_fastq'
include fasttree from "./modules/fasttree"
include filter_fastq_by_length from './modules/filter_fastq_by_length'
include mafft from './modules/mafft'
include split_reference from './modules/split_reference'
include toytree from './modules/toytree'
include snp_sites from './modules/snp_sites'

/************************** 
* SUB WORKFLOWS
**************************/

workflow artic_nCov19_wf {
    take:   
        fastq
    main:   
        artic(filter_fastq_by_length(fastq))
    emit:   
        artic.out
}


// TODO: rewrite this one to the nextstrain pipeline

workflow create_tree_snippy_wf {
    take: 
        fasta       // the nCov fasta
        references  // multiple references to compare against
    main:
        split_reference(references)

        // val(fasta), path(fasta), val(reference_name), path(references)
        snippy_input =  fasta   .combine(split_reference.out
                                .flatten()
                                .map { file -> tuple( file.baseName, file ) } 
                                )

        snippy(snippy_input)

        input_snippy_msa =  snippy.out
                                .groupTuple()
                                .map { it -> tuple(it[0], it[1][0], it[2]) }

        fasttree(
            snp_sites(
                    snippy_msa(input_snippy_msa)))
    emit:
        fasttree.out
}

// TODO: get fastaname and carry it as env / val to the toytree highlight
// this way i can highlight all samples in there
workflow create_tree_mafft_wf {
    take: 
        fasta
        references
    main:
        fasttree(
            snp_sites(
                mafft (fasta, references)))
            
    emit:
        fasttree.out
}

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
    if (params.artic_ncov19 && params.dir) { artic_nCov19_wf(cat_fastq(dir_input_ch)); fasta_input_ch = artic_nCov19_wf.out }
    if (params.artic_ncov19 && params.fastq) { artic_nCov19_wf(fastq_input_ch); fasta_input_ch = artic_nCov19_wf.out}

// analyse genome to references
    if (params.references && (params.fastq || params.fasta || params.dir)) { 
        if (!params.snippy) { create_tree_mafft_wf (fasta_input_ch, reference_input_ch); newick = create_tree_mafft_wf.out } 
        if (params.snippy)  { create_tree_snippy_wf (fasta_input_ch, reference_input_ch); newick = create_tree_snippy_wf.out }

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
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run replikation/nCov ${c_blue}--artic_ncov19${c_reset} ${c_green}--fastq 'sample_01.fasta.gz'${c_reset} \\ 
            ${c_blue}--references references.fasta${c_reset} -profile local,docker

    ${c_yellow}Workflow options:${c_reset}
    ${c_blue} --artic_ncov19 ${c_reset} one sample per fastq or fastq.gz file ${c_green}[--fastq]${c_reset} 
                                        or one dir containing fastq files of one ont run ${c_green}[--dir]${c_reset} 
    ${c_dim} Options:${c_reset} 
    ${c_dim} --primerV${c_reset}        artic-ncov2019 primer_schemes [default: ${params.primerV}]
    ${c_dim} --minLength${c_reset}      min length filter raw reads [default: ${params.minLength}]
    ${c_dim} --maxLength${c_reset}      max length filter raw reads [default: ${params.maxLength}]
 
    ${c_yellow}Analysis options:${c_reset}
    ${c_blue} --references ${c_reset}   [--references references.fasta] fasta file(s) to compare against sample 
                                        Provide sample(s) via ${c_green}[--fasta]${c_reset}
                                        Or use ${c_green}[--fastq]${c_reset} or ${c_green}[--dir]${c_reset} && ${c_blue}--<workflow>${c_reset}

    ${c_dim} Options:${c_reset} 
    ${c_dim} --snippy${c_reset}         use snippy instead of mafft to build tree

    ${c_reset}Options:
    --cores             max cores for local use [default: $params.cores]
    --memory            available memory [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    Profile:
    -profile                 local,docker -> merge profiles e.g. -profile local,docker ${c_reset}
    """.stripIndent()
}

/*

// NEXTSTRAIN container via: (?)
// https://nextstrain.org/docs/getting-started/local-installation

 https://github.com/blab/sars-like-cov

 snake file
https://github.com/blab/sars-like-cov/blob/master/Snakefile

https://nextstrain-augur.readthedocs.io/en/stable/usage/cli/tree.html
the TLDR


// AUGUR approach

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
        """
    


rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.align.output.alignment
    output:
        alignment = "results/masked.fasta"
    params:
        mask_from_beginning = 15,
        mask_from_end = 15,
        mask_sites = 18460
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """


*/