#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* Nextflow -- nCov Analysis Pipeline
* Author: christian.jena@gmail.com
*/

/************************** 
* HELP messages & checks
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
println "CPUs to use: $params.cores"
println "Memory in GB: $params.memory"
if (!params.single && params.dir) { println "Barcodes: True" }
if (params.dir || params.fastq) { println "Primerscheme: $params.primerV"  }
println "Output dir: $params.output\u001B[0m"
println " "

// profile helps
    if ( workflow.profile == 'standard' ) { exit 1, "NO EXECUTION PROFILE SELECTED, use e.g. [-profile local,docker]" }
    if (params.profile) { exit 1, "--profile is WRONG use -profile" }
    if (
        workflow.profile.contains('singularity') ||
        workflow.profile.contains('docker')
        ) { "engine selected" }
    else { println "No engine selected:  -profile EXECUTER,ENGINE" 
           println "using native installations" }
    if (
        workflow.profile.contains('local')
        ) { "executer selected" }
    else { exit 1, "No executer selected:  -profile EXECUTER,ENGINE" }

// params help
    if (!params.fasta &&  !params.dir &&  !params.fastq ) {
        exit 1, "input missing, use [--fasta] [--fastq] or [--dir]"}
    if ((params.fasta && ( params.fastq || params.dir )) || ( params.fastq && params.dir )) {
        exit 1, "To much inputs: please us either: [--fasta], [--fastq] or [--dir]"} 
    if (!params.metadata) { println "\033[0;33mNo [--metadata] file specified, skipping tree build\u001B[0m" }

/************************** 
* INPUTs
**************************/

// fasta input 
    if (params.fasta) { fasta_input_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
    }

// references input 
    if (params.references) { reference_input_ch = Channel
        .fromPath( params.references, checkIfExists: true)
    }

// metadata input 
    if (params.metadata) { metadata_input_ch = Channel
        .fromPath( params.metadata, checkIfExists: true)
    }

// fastq input
    if (params.fastq) { fastq_input_ch = Channel
        .fromPath( params.fastq, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
    }

// dir input
    if (params.dir) { dir_input_ch = Channel
        .fromPath( params.dir, checkIfExists: true, type: 'dir')
        .map { file -> tuple(file.name, file) }
    }

/************************** 
* DATABASES
**************************/

workflow build_database_wf {
    main:
        fasta_DB = Channel.fromPath( workflow.projectDir + "/database/ena_*.fasta" , checkIfExists: true)
        text_DB = Channel.fromPath( workflow.projectDir + "/database/ena_*.txt", checkIfExists: true)
    
        create_database(fasta_DB, text_DB)
    emit:
        create_database.out[0]
        create_database.out[1]
}

/************************** 
* MODULES
**************************/

include artic from './modules/artic' 
include augur_align from './modules/augur'
include augur_tree from './modules/augur'
include augur_tree_refine from './modules/augur'
include bwa_samtools from './modules/bwa_samtools'
include coverage_plot from './modules/coverage_plot'
include create_database from './modules/create_database'
include filter_fastq_by_length from './modules/filter_fastq_by_length'
include guppy_gpu from './modules/guppy'
include mask_alignment from './modules/mask_alignment'
include pangolin from './modules/pangolin' 
include quality_genome_filter from './modules/quality_genome_filter'
include toytree from './modules/toytree'

/************************** 
* SUB WORKFLOWS
**************************/

workflow basecalling_wf {
    take: 
        dir_input  
    main:
        
        guppy_gpu(dir_input)
        
        if (params.single) { fastq_channel = guppy_gpu.out }

        else { fastq_channel = guppy_gpu.out
                            .map { it -> it[1] }
                            .flatten()
                            .map { it -> [ it.simpleName, it ] }
            }
    
    emit:
        fastq_channel
} 


workflow artic_nCov19_wf {
    take:   
        fastq
    main: 

        // assembly  
        artic(filter_fastq_by_length(fastq))

        // validate fasta
        coverage_plot(
            bwa_samtools(
                artic.out.join(filter_fastq_by_length.out)))

    emit:   
        artic.out
}


workflow create_tree_wf {
    take: 
        fasta       // the nCov fasta (own samples or reconstructed here)
        references  // multiple references to compare against
        metadata    // tsv file of meta data  strain country date
    main:

        align_reference = Channel.fromPath( workflow.projectDir + "/data/reference_nCov19/MN908947.gb", checkIfExists: true)

        quality_genome_filter(fasta)

        // from [val, file] to [files]
        collect_fasta = quality_genome_filter.out[0].map{ it -> it [1]}
                                                    .collect()

        augur_tree(
            mask_alignment(
                augur_align(collect_fasta, references, align_reference)))

        augur_tree_refine(augur_tree.out, metadata)

    emit:
        augur_tree_refine.out
}

/*
 TODO: get fastaname and carry it as env / val to the toytree highlight
 this way i can highlight all samples in there
 could be maybe done with a channel extracting the names which merges in here

 also highlight location data or annotate this to the nodes
*/

workflow toytree_wf {
    take: 
        trees  
    main:
        toytree(trees)
    emit:
        toytree.out
} 

workflow determine_lineage_wf {
    take: 
        fasta  
    main:
        pangolin(fasta)

        // collect lineage also to a summary     
        channel_tmp = pangolin.out[1]
                .splitCsv(header: true, sep: ',')
                .collectFile(seed: 'taxon,lineage,SH-alrt,UFbootstrap,lineages_version,status,note\n', 
                            storeDir: params.output + "/summary/") {
                            row -> [ "metadata.tsv", row.taxon + ',' + row.lineage + ',' + row.'SH-alrt' + ',' + 
                            row.'UFbootstrap' + ',' + row.lineages_version + ',' + row.status + ',' + row.note + '\n']
                            }

    emit:
        pangolin.out[0]
} 

/************************** 
* MAIN WORKFLOW
**************************/

workflow {

// 1. reconstruct genomes
    if (params.dir) { artic_nCov19_wf(basecalling_wf(dir_input_ch)); fasta_input_ch = artic_nCov19_wf.out }
    if (params.fastq) { artic_nCov19_wf(fastq_input_ch); fasta_input_ch = artic_nCov19_wf.out}

// 2. analyse genomes to references and build tree
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
    }

// 3. plot tree
    if (params.metadata) { toytree_wf(newick) }

// 4. determine lineage
    if (params.fastq || params.fasta || params.dir) {
        determine_lineage_wf(fasta_input_ch)

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
    
    ${c_green}poreCov${c_reset} | A Nextflow nCov19 workflow for nanopore data
    
    ${c_yellow}Usage examples:${c_reset}
    nextflow run replikation/poreCov --fastq 'sample_01.fasta.gz' --cores 14 -profile local,singularity

    ${c_yellow}Inputs (choose one):${c_reset}
    --dir           one fast5 dir of a nanopore run containing multiple samples (barcoded)
                    [--dir fast5/ --single] if you only one sample (no barcodes)
                    ${c_dim}(not implemented yet) autorename barcodes via [--barcodeIDs rename.csv] 
                        Per line: 01;samplename${c_reset}
                    ${c_dim}[basecalling - demultiplexing - nCov genome reconstruction]${c_reset}

    --fastq         one fastq or fastq.gz file per sample or
                    multiple file-samples: --fastq 'sample_*.fasta.gz'
                    ${c_dim}[nCov genome reconstruction]${c_reset}

    --fasta         direct input of genomes, one file per genome
                    ${c_dim}(multifasta not implemented yet) --multifasta ${c_reset}
                    ${c_dim}[phylogenetic tree contruction]${c_reset}

    ${c_yellow}Parameters - Basecalling${c_reset}
    --localguppy    use a native guppy installation instead of a gpu-guppy-docker 
                    native guppy installation is used by default for singularity or conda

    ${c_yellow}Parameters - nCov genome reconstruction${c_reset}
    --primerV       artic-ncov2019 primer_schemes [default: ${params.primerV}]
    --minLength     min length filter raw reads [default: ${params.minLength}]
    --maxLength     max length filter raw reads [default: ${params.maxLength}]

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

    ${c_yellow}Options:${c_reset}
    --cores         max cores for local use [default: $params.cores]
    --memory        available memory [default: $params.memory]
    --output        name of the result folder [default: $params.output]
    --cachedir      defines the path where singularity images are cached
                    [default: $params.cachedir] 

    ${c_yellow}Execution/Engine profiles:${c_reset}
    nCov supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.:
     -profile ${c_green}local${c_reset},${c_blue}docker${c_reset}

      ${c_green}Executer${c_reset} (choose one):
      local
      ${c_blue}Engines${c_reset} (choose one):
      docker
      singularity
      conda (not yet implemented )
    """.stripIndent()
}
