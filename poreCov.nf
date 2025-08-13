#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* Nextflow -- SARS-CoV-2 Analysis Pipeline
* Author: christian.jena@gmail.com
*/

/************************** 
* MODULES
**************************/

include { get_fast5 } from './modules/get_fast5_test_data.nf'
include { get_nanopore_fastq } from './modules/get_fastq_test_data.nf'
include { get_fasta } from './modules/get_fasta_test_data.nf'
include { align_to_reference } from './modules/align_to_reference.nf'
include { split_fasta } from './modules/split_fasta.nf'
include { filter_fastq_by_length } from './modules/filter_fastq_by_length.nf'
include { count_mixed_sites } from './modules/count_mixed_sites.nf'

/************************** 
* Workflows
**************************/

include { artic_ncov_wf } from './workflows/artic_nanopore_nCov19.nf'
include { basecalling_wf } from './workflows/basecalling.nf'
include { collect_fastq_wf } from './workflows/collect_fastq.nf'
include { create_json_entries_wf } from './workflows/create_json_entries.nf'
include { create_summary_report_wf } from './workflows/create_summary_report.nf'
include { determine_lineage_wf } from './workflows/determine_lineage.nf'
include { determine_mutations_wf } from './workflows/determine_mutations.nf'
include { genome_quality_wf } from './workflows/genome_quality.nf'
include { read_classification_wf; read_screening_freyja_wf; read_screening_lsc_wf} from './workflows/read_classification'
include { read_qc_wf } from './workflows/read_qc.nf'
include { rki_report_wf } from './workflows/provide_rki.nf'
include { pangolin } from './workflows/process/pangolin.nf'

/************************** 
* Begin of Workflow
**************************/
workflow {

    header()

/************************** 
* HELP messages & checks
**************************/

// try to check for poreCov releases

    println " "
    println "  Latest available poreCov release: " + Git_PoreCov.getLatestVersion()
    println "  If neccessary update via: nextflow pull replikation/poreCov"
    println "________________________________________________________________________________"


// Log infos based on user inputs
    if ( params.help ) { exit 0, helpMSG() }

// profile helps
    if ( workflow.profile == 'standard' ) { exit 1, "NO EXECUTION PROFILE SELECTED, use e.g. [-profile local,docker]" }
    if (params.profile) { exit 1, "--profile is WRONG use -profile" }
    if (
        workflow.profile.contains('singularity') ||
        workflow.profile.contains('nanozoo') ||
        workflow.profile.contains('ukj_cloud') ||
        workflow.profile.contains('stub') ||
        workflow.profile.contains('docker')
        ) { "engine selected" }
    else { println "No engine selected:  -profile EXECUTER,ENGINE" 
           println "using native installations" }
    if (
        workflow.profile.contains('nanozoo') ||
        workflow.profile.contains('ukj_cloud') ||
        workflow.profile.contains('local') ||
        workflow.profile.contains('stub') ||
        workflow.profile.contains('slurm')
        ) { "executer selected" }
    else { exit 1, "No executer selected:  -profile EXECUTER,ENGINE" }

    if (workflow.profile.contains('local')) {
        println "\033[2m Using $params.cores/$params.max_cores CPU threads per process for a local run [--max_cores]\u001B[0m"
        println " "
    }
    if ( workflow.profile.contains('singularity') ) {
        println ""
        println "\033[0;33mWARNING: Singularity image building sometimes fails!"
        println "Multiple resumes (-resume) and --max_cores 1 --cores 1 for local execution might help.\033[0m\n"
    }

// params help
    if (!workflow.profile.contains('test_fastq') && !workflow.profile.contains('test_fast5') && !workflow.profile.contains('test_fasta')) {
        if (!params.fasta &&  !params.fast5 &&  !params.fastq && !params.fastq_pass ) {
            exit 1, "input missing, use [--fasta] [--fastq] [--fastq_pass] or [--fast5]" }
        if ( params.fastq && params.fastq_pass ) { exit 1, "Please use either: [--fastq] or [--fastq_pass]"}
        if ( params.fasta && ( params.fastq || params.fast5 || params.fastq_pass)) { exit 1, "Please use [--fasta] without inputs like: [--fastq], [--fastq_pass], [--fast5]" }
        if (params.list && params.fasta) { exit 1, "[--fasta] and [--list] is not supported" }


    }
    if ( (params.cores.toInteger() > params.max_cores.toInteger()) && workflow.profile.contains('local')) {
            exit 1, "More cores (--cores $params.cores) specified than available (--max_cores $params.max_cores)" }

    if ( params.single && params.samples ) { exit 1, "Sample input [--samples] not supported for [--single]" }

// check that input params are used as such
    if (params.fasta == true) { exit 5, "Please provide a fasta file via [--fasta]" }
    if (params.fastq == true) { exit 5, "Please provide a fastq files (one per sample) via [--fastq]" }
    if (params.fastq_pass == true) { exit 5, "Please provide a fastq_pass dir via [--fastq_pass]" }
    if (params.fast5 == true) { exit 5, "Please provide a fast5 dir via [--fast5]" }
    if (params.minLength && !params.minLength.toString().matches("[0-9]+")) { exit 5, "Please provide an integer number (e.g. 300) as minimal read length via [--minLength]" }
    if (params.maxLength && !params.maxLength.toString().matches("[0-9]+")) { exit 5, "Please provide an integer number (e.g. 300) as maximum read length via [--maxLength]" }
    if (params.nanopolish)   { println "\033[0;33mWarning: Parameter [--nanopolish] is deprecated, ignoring flag.\033[0m" }
    if (params.medaka_model) { println "\033[0;33mWarning: Parameter [--medaka_model] is deprecated, please use [--clair3_model_dir] and / or [--clair3_model_name] to specify a non default model.\033[0m" }


// check correct usage of param-flags
    if (params.extended && !params.samples ) { exit 5, "When using --extended you need to specify also a sample.csv via [--samples]" }
    if (!params.freyja == true && !params.freyja == false) {exit 5, "Please provide no input to [--freyja]"}
    if (!params.lcs == true && !params.lcs == false) {exit 5, "Please provide no input to [--lcs]"}
    if (params.screen_reads && !params.lcs && !params.freyja) {exit 5, "When using [--screen_reads] you also need to use at least one: [--freyja] or [--lcs]"}
    if (!params.screen_reads && params.lcs) {exit 5, "[--lcs] requires [--screen_reads] to work"}
    if (!params.screen_reads && params.freyja) {exit 5, "[--freyja] requires [--screen_reads] to work"}


        
// validate primer scheme version format
    def fetched_version = "${params.primerV}" =~ /(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)/
    def legacy_primerV = ['V1','V1200','V2','V3','V4','V4.1','V5','V5.1','V5.2.0_1200','V5.3.2_400']
    if (!fetched_version && !("${params.primerV}".contains('.bed')) && !(legacy_primerV.any{params.primerV.contains(it)})){ exit 1, "Invalid scheme version format '${params.primerV}' provided, please provide a version in the format 'vX.X.X', e.g. v1.0.0" }
    if ("${params.primerV}".contains('.bed') && !params.primerRef){ exit 1, "Custom primer scheme '${params.primerV}' was provided without primer reference. Please pass a primer scheme reference sequence via [--primerRef]." }
    
    if ("${params.primerV}".contains('_') && !("${params.primerV}".contains('.bed')) || (legacy_primerV.any{params.primerV.contains(it)} && !("${params.primerV}".contains('.bed')))){ 
        println "\033[2mPrimer scheme version in legacy format detected. Using local nCov-19 primer schemes.\033[0m"
        legacy_primerV = true 
    } else {
        legacy_primerV = false
    }

// validating sample table
    if (params.samples) {  

        // check that the rows _id and Status can be found
        // checks afterwards that no fields are empty
        Channel.fromPath( params.samples, checkIfExists: true)
            .splitCsv(header: false, sep: ',')
            .take( 1 )
            .map { row ->  
                if ( !("_id" in row) ) { exit 6, "The column '_id' was not found in $params.samples, hidden symbols? Use a editor to generate the csv file" }
                if ( !("Status" in row) ) { exit 6, "The column 'Status' was not found in $params.samples" }
            }
            .mix(
            Channel.fromPath( params.samples, checkIfExists: true)
                .splitCsv(header: true, sep: ',')
                .map { row -> 
                    if (!row.'Status') { exit 6, "A Status field appears to be empty in the file $params.samples" }
                    if (!row.'_id') { exit 6, "A _id field appears to be empty in the file $params.samples"} 
                }
            )
    }


/************************** 
* INPUTs
**************************/

// fasta input 
    if (params.fasta && !workflow.profile.contains('test_fasta')) { 
        fasta_input_raw_ch = Channel
        .fromPath( params.fasta, checkIfExists: true)
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

// fastq input or via csv file
    if (params.fastq.toString().contains('.csv') && !workflow.profile.contains('test_fastq')) { 
        fastq_file_ch = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
    }
    else if (params.fastq && !workflow.profile.contains('test_fastq')) { 
        fastq_file_ch = Channel
        .fromPath( params.fastq, checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
    }

// fastq raw input direct from basecalling
    if (params.fastq_pass && !workflow.profile.contains('test_fastq')) { 
        fastq_dir_ch = Channel
        .fromPath( params.fastq_pass, checkIfExists: true, type: 'dir')
        .map { file -> tuple(file.simpleName, file) }
    }

// dir input
    if (params.fast5 && !workflow.profile.contains('test_fast5')) { dir_input_ch = Channel
        .fromPath( params.fast5, checkIfExists: true, type: 'dir')
        .map { file -> tuple(file.name, file) }
    }

// samples input 
    if (params.samples) { 
        samples_input_ch = Channel.fromPath( params.samples, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .map { row -> tuple ("barcode${row.Status[-2..-1]}", "${row._id.replace( " ", "")}")}

        samples_file_ch = Channel.fromPath( params.samples, checkIfExists: true)
            .splitCsv(header: true, sep: ',')
            .collectFile(seed: '_id,Status\n') {
                        row -> [ "input.csv", row.'_id'.replace( " ", "") + ',' + row.'Status'.replace( " ", "") + '\n']
                        }
    }

    else { samples_file_ch = Channel.from( ['deactivated'] ) }

    // extended input
    if (params.samples && params.extended) { 
        extended_input_ch = Channel.fromPath( params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .collectFile() {
                    row -> [ "extended.csv", row.'_id'.replace( " ", "") + ',' + row.'Submitting_Lab' + ',' + row.'Isolation_Date' + ',' + 
                    row.'Seq_Reason' + ',' + row.'Sample_Type'.replace( " ", "") + '\n']
                    }
    }
    
    else { extended_input_ch = Channel.from( ['deactivated', 'deactivated'] ) }

/************************** 
* Automatic Pangolin version updates, with fail save
**************************/

    if (params.update) {
    println "\033[0;33mWarning: Running --update might not be poreCov compatible!\033[0m"
        if ( Dockerhub_Pangolin.IsPangoAvailable().toString() == "true" ) { 
            tagname = 'https://registry.hub.docker.com/v2/repositories/nanozoo/pangolin-v4/tags/'.toURL().text.split(',"name":"')[1].split('","')[0]
            pangolindocker = "nanozoo/pangolin-v4:" + tagname
            println "\033[0;32mFound latest pangolin container, using: " + pangolindocker + " \033[0m" 

            tagname = 'https://registry.hub.docker.com/v2/repositories/nanozoo/nextclade3/tags/'.toURL().text.split(',"name":"')[1].split('","')[0]
            nextcladedocker = "nanozoo/nextclade3:" + tagname 
            println "\033[0;32mFound latest nextclade3 container, using: " + nextcladedocker + " \033[0m"
        } 
        else { 
            println "\033[0;33mCould not find the latest pangolin container, trying: " + params.defaultpangolin + "\033[0m"
            pangolindocker = params.defaultpangolin 

            println "\033[0;33mCould not find the latest nextclade3 container, trying: " + params.defaultnextclade + "\033[0m"
            nextcladedocker = params.defaultnextclade 
        } 
    }
    else { pangolindocker = params.defaultpangolin ; nextcladedocker = params.defaultnextclade  }

    if ( params.screen_reads && params.lcs_ucsc_update ){
        if ( Hgdownload_lcs.IsLcsAvailable().toString() == "true" ) { 
            lsc_ucsc_work_version = 'https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt'.toURL().text.split('\\(')[1].split('\\)')[0]
            log.info "\033[0;32mFound latest UCSC version, using: " + lsc_ucsc_work_version + " \033[0m" 
        }
        else { 
            log.info "\033[0;33mCould not find the latest UCSC version, trying: " + params.lcs_ucsc_version + "\033[0m"
            lsc_ucsc_work_version = params.lcs_ucsc_version
        }
    } else { lsc_ucsc_work_version = params.lcs_ucsc_version}

    defaultMSG()
    if ( params.fast5 || workflow.profile.contains('test_fast5') ) { basecalling() }
    if (!params.fasta && !workflow.profile.contains('test_fasta')) { read_length() }
    rki()

/************************** 
* MAIN WORKFLOW
**************************/

    // 0. Test profile data
        if ( workflow.profile.contains('test_fast5')) { dir_input_ch =  get_fast5().map {it -> ['SARSCoV2', it] } }
        if ( workflow.profile.contains('test_fastq')) { fastq_input_raw_ch =  get_nanopore_fastq().map {it -> ['SARSCoV2', it] } }
        if ( workflow.profile.contains('test_fasta')) { fasta_input_raw_ch =  get_fasta() }

    // 1. Reconstruct genomes
        // fast5
        if ( (params.fast5 && !params.fastq && !params.fastq_pass) || workflow.profile.contains('test_fast5')) {
            basecalling_wf(dir_input_ch)
            
            // rename barcodes
                if (params.samples) { 
                    fastq_from5_ch = basecalling_wf.out[0].join(samples_input_ch).map { it -> tuple(it[2],it[1]) }
                    basecalling_wf.out[0].join(samples_input_ch).ifEmpty{ exit 2, "Could not match barcode numbers from $params.samples to the read files, some typo?"} 
                    }
                else if (!params.samples) { fastq_from5_ch = basecalling_wf.out[0] }
            
            filtered_reads_ch = filter_fastq_by_length(fastq_from5_ch)
            noreadsatall = filtered_reads_ch.ifEmpty{ log.info "\033[0;33mNot enough reads in all samples, please investigate $params.output/$params.readqcdir\033[0m" }
            read_classification_wf(filtered_reads_ch)

            artic_ncov_wf(legacy_primerV, filtered_reads_ch, params.artic_normalize) 
            fasta_input_ch = artic_ncov_wf.out.assembly

            // count mixed sites
            if (params.primerV.toString().contains(".bed")) {
                external_primer_schemes = artic_ncov_wf.out.primer_dir
            }
            else {
                external_primer_schemes = file(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            }
            count_mixed_sites(artic_ncov_wf.out.trimmed_bam.join(artic_ncov_wf.out.vcf).join(artic_ncov_wf.out.failed_vcf), external_primer_schemes)
        }
        // fastq input via dir and or files
        if ( (params.fastq || params.fastq_pass) || workflow.profile.contains('test_fastq')) { 
            if (params.fastq_pass && !params.fastq) { fastq_input_raw_ch = collect_fastq_wf(fastq_dir_ch) }
            if (!params.fastq_pass && params.fastq) { fastq_input_raw_ch = fastq_file_ch }

            // rename barcodes based on --samples input.csv
                if (params.samples) { fastq_input_ch = fastq_input_raw_ch.join(samples_input_ch).map { it -> tuple(it[2],it[1])} 
                   fastq_input_raw_ch.join(samples_input_ch).ifEmpty{ exit 2, "Could not match barcode numbers from $params.samples to the read files, some typo?"} 
                }
                else if (!params.samples) { fastq_input_ch = fastq_input_raw_ch }

            read_qc_wf(fastq_input_ch)
            filtered_reads_ch = filter_fastq_by_length(fastq_input_ch)
            noreadsatall = filtered_reads_ch.ifEmpty{ log.info "\033[0;33mNot enough reads in all samples, please investigate $params.output/$params.readqcdir\033[0m" }
            read_classification_wf(filtered_reads_ch)

            // genome reconstruction with artic
            artic_ncov_wf(legacy_primerV, filtered_reads_ch, params.artic_normalize)
            fasta_input_ch = artic_ncov_wf.out.assembly
            
            // count mixed_sites
            if (params.primerV.toString().contains(".bed")) {
                external_primer_schemes = artic_ncov_wf.out.primer_dir
            }
            else {
                external_primer_schemes = file(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            }
            count_mixed_sites(artic_ncov_wf.out.trimmed_bam.join(artic_ncov_wf.out.vcf).join(artic_ncov_wf.out.failed_vcf), external_primer_schemes)
        }

    // 2. Genome quality, lineages, clades and mutations
        // fasta input
        if ( params.fasta || workflow.profile.contains('test_fasta' ) ) {
            fasta_input_ch = split_fasta(fasta_input_raw_ch).flatten().map { it -> tuple(it.simpleName, it) }
        }

        determine_lineage_wf(fasta_input_ch, pangolindocker)
        determine_mutations_wf(fasta_input_ch, nextcladedocker)
        genome_quality_wf(fasta_input_ch, reference_for_qc_input_ch)

    // 3. Specialised outputs (rki, json)
        rki_report_wf(genome_quality_wf.out.president_valid, genome_quality_wf.out.president_invalid, extended_input_ch)

        if (params.samples) {
            create_json_entries_wf(determine_lineage_wf.out, genome_quality_wf.out.president_valid, determine_mutations_wf.out)
        }

    // 4. Summary output
        if (params.fasta || workflow.profile.contains('test_fasta')) {
            taxonomic_read_classification_ch = Channel.from( ['deactivated', 'deactivated', 'deactivated'] ).collect()
            alignments_ch = Channel.from( ['deactivated'] )

        } else {
            taxonomic_read_classification_ch = read_classification_wf.out.kraken
            if (params.screen_reads) {
                if (params.lcs) {
                    read_screening_lsc_wf(filtered_reads_ch, lsc_ucsc_work_version)
                }
                if (params.freyja) {
                    read_screening_freyja_wf(artic_ncov_wf.out.trimmed_bam.map{it -> [it[0], it[1]]}.combine(reference_for_qc_input_ch))
                }
            }
            alignments_ch = align_to_reference(filtered_reads_ch.combine(reference_for_qc_input_ch))
        }
        if (params.fasta || workflow.profile.contains('test_fasta')) {
            alt_allele_ratio_ch = Channel.from( ['deactivated'] )
        } else {
            alt_allele_ratio_ch = count_mixed_sites.out.stats
        }

/*
        if (params.samples) {
            samples_table_ch = Channel.fromPath( params.samples, checkIfExists: true)
        }
        else { samples_table_ch = Channel.from( ['deactivated'] ) }
*/
        create_summary_report_wf(determine_lineage_wf.out, genome_quality_wf.out.president_valid, determine_mutations_wf.out,
                                taxonomic_read_classification_ch, alt_allele_ratio_ch, alignments_ch, samples_file_ch, nextcladedocker)

}

/*************  
* --help
*************/
def helpMSG() {
    def c_green = "\033[0;32m";
    def c_reset = "\033[0m";
    def c_yellow = "\033[0;33m";
    def c_blue = "\033[0;34m";
    def c_dim = "\033[2m";
    log.info """
 ${c_dim}  ${c_reset}    
\033[0;33mUsage examples:${c_reset}

    nextflow run replikation/poreCov --update --fastq '*.fasta.gz' -r 1.3.0 -profile local,singularity

${c_yellow}Inputs (choose one):${c_reset}
  --fast5         One fast5 dir of a nanopore run containing multiple samples (barcoded)
                  Add the flag [--single] if no barcodes were used
                  ${c_dim}[Basecalling + Genome reconstruction + Lineage + Reports]${c_reset}

  --fastq         One fastq or fastq.gz file per sample or
                  Multiple file-samples: --fastq 'sample_*.fastq.gz'.
                  Or Table .csv input: ${c_dim}Sample_name,Absolut_path${c_reset} per line, NO header or SPACE
                  ${c_dim}[Genome reconstruction + Lineage + Reports]${c_reset}

  --fastq_pass    The fastq_pass dir after bascalling with barcodes
                  --fastq_pass 'fastq_pass/'
                  Add the flag [--single] if no barcodes were used
                  ${c_dim}[Genome reconstruction + Lineage + Reports]${c_reset}

  --fasta         Direct input of genomes - supports multi-fasta file(s) - can be gzip compressed (.gz)
                  ${c_dim}[Lineage + Reports]${c_reset}

${c_yellow}Workflow control (optional)${c_reset}
  --artic_normalize   Normalise reads to moderate coverage, saves runtime [default: $params.artic_normalize]
                      ${c_dim}(after mapping and before variant calling in the ARTIC bioinformatics pipeline)
                      Use `--artic_normalize False` to turn off this normalisation.${c_reset}
  --update            Always try to use latest pangolin & nextclade release [default: $params.update]
  --samples           Table .csv input (header: Status,_id), renames barcodes (Status) by name (_id), e.g.:
                          ${c_dim}Status,_id
                          barcode01,sample2011XY
                          BC02,thirdsample_run${c_reset}
  --extended          poreCov utilizes from --samples these additional headers:
                          ${c_dim}Submitting_Lab,Isolation_Date,Seq_Reason,Sample_Type${c_reset}
  --nanopolish        Use nanopolish instead of medaka for ARTIC (needs --fast5)
                          to skip basecalling provide --fastq or --fastq_pass and a sequencing_summary.txt
                          e.g: --fastq_pass fastq_pass/ --nanopolish sequencing_summary.txt --fast5 fast5_dir/
  --screen_reads      Determines Pangolineage of each read (takes time, needs --freyja and/or --lcs)
  --scorpio           Skip Scorpio in pangolin run [default: $params.scorpio]
                          ${c_dim}Scorpio overwrites Usher results, might lead to many unassigned samples${c_reset}

${c_yellow}Parameters - Lineage detection on reads (see --screen_reads, optional)${c_reset}
  --freyja            Activates --screen_reads via freyja
  --freyja_update     Update freyja's barcode-db prior to running
  
  --lcs               Activate --screen_reads via lcs
  --lcs_ucsc_version  Create marker table based on a specific UCSC SARS-CoV-2 tree (e.g. '2022-05-01'). 
                      Predefined (default) uses the repository table (not up-to-date) [default: $params.lcs_ucsc_version]
                      ${c_dim}Available trees: https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2${c_reset}
  --lcs_cutoff        Plot linages above this threshold [default: $params.lcs_cutoff]       
  --lcs_ucsc_update   Use latest UCSC SARS-CoV-2 tree. Overwrites --lcs_ucsc_version [default: $params.lcs_ucsc_update]
        ${c_dim}From: https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt${c_reset}
  --lcs_ucsc_downsampling  Downsample sequences for --lcs_ucsc_update. Use 'None' to turn off [default: $params.lcs_ucsc_downsampling]
                           ${c_dim}Attention! High Resource and MEMORY usage without downsampling${c_reset}
  --lcs_variant_groups  Variant groups table (.tsv) for --lcs_ucsc_update.
                        Use 'default' for predefined groups. [default: $params.lcs_variant_groups]
        ${c_dim}From: https://github.com/rki-mf1/LCS/blob/master/data/variant_groups.tsv${c_reset} 

${c_yellow}Parameters - Basecalling  (optional)${c_reset}
  --localguppy    Use a native installation of guppy instead of a gpu-docker or gpu_singularity 
  --guppy_cpu     Use cpus instead of gpus for basecalling ${c_dim}(super slow)${c_reset}
  --one_end       Removes the recommended "--require_barcodes_both_ends" from guppy demultiplexing
                  ${c_dim}try this if to many barcodes are unclassified (beware - results might not be trustworthy)${c_reset}
  --guppy_model   Guppy basecalling model [default: ${params.guppy_model}]
                  ${c_dim}e.g. "dna_r9.4.1_450bps_hac.cfg" or "dna_r9.4.1_450bps_sup.cfg"${c_reset}

${c_yellow}Parameters - SARS-CoV-2 genome reconstruction (optional)${c_reset}
  --primerV             Supported primer variants or primer bed files - choose one [default: ${params.primerV}]
                            ${c_dim}ARTIC (>v1.6.0) :${c_reset} V1, V2, V3, V4, V4.1, V.5, V.5.1, V.5.3.2_400
                            ${c_dim}ARTIC (>=v1.6.0):${c_reset} v1.0.0, v2.0.0, v3.0.0, v4.0.0, v4.1.0, v5.0.0, v5.1.0, v5.3.2
                            ${c_dim}NEB:${c_reset} VarSkipV1a, VarSkipV2, VarSkipV2b
                            ${c_dim}Other:${c_reset} V1200, V5.2.0_1200 ${c_dim}(also known as midnight)${c_reset}
                            ${c_dim}Primer bed file:${c_reset} e.g. primers.bed  ${c_dim}See Readme for more help${c_reset}
  --schemeLength          primer scheme length, e.g. 400, 700; artic remote primers are length 400, varvamp remote primers 700 [default: ${params.schemeLength}]
  --rapid                 rapid-barcoding-kit was used [default: ${params.rapid}]
  --minLength             min length filter raw reads [default: 100]
  --maxLength             max length filter raw reads [default: 700 (primer-scheme: V1-4, rapid); 1500 (primer-scheme: V1200, V5.2.0_1200)]
  --min_depth             nucleotides below min depth will be masked to "N" [default ${params.min_depth}]
  --clair3_model_dir      directory to look for clair3 model files [default: ${params.clair3_model_dir}]
  --clair3_model_name     clair3 model for the artic workflow [default: ${params.clair3_model_name}]

${c_yellow}Parameters - Genome quality control  (optional)${c_reset}
  --reference_for_qc      Reference FASTA for consensus qc (optional, wuhan is provided by default)
  --seq_threshold         Global pairwise ACGT sequence identity threshold [default: ${params.seq_threshold}] 
  --n_threshold           Consensus sequence N threshold [default: ${params.n_threshold}] 

${c_yellow}Options  (optional)${c_reset}
  --cores         Amount of cores for a process (local use) [default: $params.cores]
  --max_cores     Max amount of cores for poreCov to use (local use) [default: $params.max_cores]
  --memory        Available memory [default: $params.memory]
  --output        Name of the result folder [default: $params.output]
  --cachedir      Defines the path where singularity images are cached
                  [default: $params.cachedir]
  --krakendb      Provide a .tar.gz kraken database [default: auto downloads one]

${c_yellow}Execution/Engine profiles (choose executer and engine${c_reset}
    poreCov supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} 
    examples:
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

    ${c_dim}Note: 'singularity' automatically passes the following environment variables to the container: 
    HTTPS_PROXY, HTTP_PROXY, http_proxy, https_proxy, FTP_PROXY, ftp_proxy${c_reset}
    """.stripIndent()
}

def header(){
    def c_green = "\033[0;32m";
    def c_reset = "\033[0m";
    log.info """
________________________________________________________________________________
    
${c_green}poreCov${c_reset} | A Nextflow SARS-CoV-2 workflow for nanopore data
    """
}

def defaultMSG(){
    def c_green = "\033[0;32m";
    def c_reset = "\033[0m";
    def c_dim = "\033[2m";
    log.info """
    ${c_green}
    Profile:             $workflow.profile${c_reset}
    ${c_dim}Current User:        $workflow.userName
    Nextflow-version:    $nextflow.version
    poreCov-version:     $workflow.revision
    ${c_reset}
    Pathing:
    ${c_dim}Workdir location [-work-Dir]:
        $workflow.workDir
    Output dir [--output]: 
        $params.output
    Databases location [--databases]:
        $params.databases
    Singularity cache dir [--cachedir]: 
        $params.cachedir ${c_reset}
    
    Parameters:
    ${c_dim}Clair3 model:         $params.clair3_model_name [--clair3_model_name]
    Clair3 model dir:     $params.clair3_model_dir
    Min depth nucleotide: $params.min_depth [--min_depth]
    Latest Pangolin/Nextclade?: $params.update [--update]
    CPUs to use:          $params.cores [--cores]
    Memory in GB:         $params.memory [--memory] ${c_reset} """.stripIndent()
}

def basecalling() {
    def c_reset = "\033[0m";
    def c_dim = "\033[2m";
    log.info """
    Basecalling options:
    ${c_dim}Use local guppy?        $params.localguppy [--localguppy]  
    One end demultiplexing? $params.one_end [--one_end]
    Basecalling via CPUs?   $params.guppy_cpu [--guppy_cpu]
    Basecalling modell:     $params.guppy_model [--guppy_model]
    Rapid-barcode-kit?:     $params.rapid [--rapid] ${c_reset} """.stripIndent()
}

def rki() {
    def c_reset = "\033[0m";
    def c_dim = "\033[2m";    
    log.info """
    RKI output for german DESH upload:
    ${c_dim}Output stored at:      $params.output/$params.rkidir  
    Min Identity to NC_045512.2: $params.seq_threshold [--seq_threshold]
    Min Depth used:        $params.min_depth [--min_depth]
       Min Depth should be 20 or more for RKI upload
    Proportion cutoff N:   $params.n_threshold [--n_threshold] ${c_reset} 
    """.stripIndent()
}

def read_length() {
    def c_reset = "\033[0m";
    def c_dim = "\033[2m";  
    def log_msg_read_min_length = params.minLength
    def log_msg_read_max_length = params.maxLength

    if ( params.primerV.matches('V1200') || params.primerV.matches('V5.2.0_1200') || params.schemeLength == 1200) {

        if ( params.primerV.matches('V1200') || params.primerV.matches('V5.2.0_1200')){
            println "\033[0;33mWarning: Definition of primer scheme length via --primerV is deprecated, please use --schemeLength instead. Setting length to 1200 ... ${c_reset}"
            params.schemeLength = 1200
        }

        if ( !params.minLength ) { log_msg_read_min_length = 400 }
        if ( !params.maxLength ) { log_msg_read_max_length = 1500 }
    }
    else {
        if ( !params.minLength ) { log_msg_read_min_length = 200 }
        if ( !params.maxLength ) { log_msg_read_max_length = 700 }
    }
    if (log_msg_read_max_length < log_msg_read_min_length) {exit 5, "--maxLength ${log_msg_read_max_length} needs to be greater than --minlength ${log_msg_read_min_length}."}

    log.info """
    Primerscheme:        $params.primerV [--primerV]
    Length:              $params.schemeLength [--schemeLength]
    ${c_dim}Min read-length set to: $log_msg_read_min_length [--minLength]
    Max read-length set to: $log_msg_read_max_length [--maxLength]${c_reset}""".stripIndent()
}
