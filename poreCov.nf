#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
* Nextflow -- SARS-CoV-2 Analysis Pipeline
* Author: christian.jena@gmail.com
*/

/************************** 
* HELP messages & checks
**************************/

header()

/* 
Nextflow version check  
Format is this: XX.YY.ZZ  (e.g. 20.07.1)
change below
*/

XX = "21"
YY = "04"
ZZ = "0"

if ( nextflow.version.toString().tokenize('.')[0].toInteger() < XX.toInteger() ) {
println "\033[0;33mporeCov requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}
else if ( nextflow.version.toString().tokenize('.')[1].toInteger() == XX.toInteger() && nextflow.version.toString().tokenize('.')[1].toInteger() < YY.toInteger() ) {
println "\033[0;33mporeCov requires at least Nextflow version " + XX + "." + YY + "." + ZZ + " -- You are using version $nextflow.version\u001B[0m"
exit 1
}

/* 
try to check for poreCov releases
*/

static boolean netIsAvailable() {
    try {
        final URL url = new URL("https://api.github.com/repos/replikation/poreCov/releases/latest");
        final URLConnection conn = url.openConnection();
        conn.connect();
        conn.getInputStream().close();
        return true;
    } catch (MalformedURLException e) {
        return false;
    } catch (IOException e) {
        return false;
    }
}

def gitcheck = netIsAvailable()

if ( gitcheck.toString() == "true" ) { porecovrelease = 'https://api.github.com/repos/replikation/poreCov/releases/latest'.toURL().text.split('"tag_name":"')[1].split('","')[0] } 
if ( gitcheck.toString() == "false" ) { porecovrelease = 'Could not get version info' } 


println " "
println "  Latest available poreCov release: " + porecovrelease
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
    if (params.fastq && params.list && !workflow.profile.contains('test_fastq')) { 
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
    if (params.fastq_pass && params.list && !workflow.profile.contains('test_fastq')) { 
        fastq_dir_ch = Channel
        .fromPath( params.fastq_pass, checkIfExists: true )
        .splitCsv()
        .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true, type: 'dir')] }
    }
    else if (params.fastq_pass && !workflow.profile.contains('test_fastq')) { 
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

static boolean DockernetIsAvailable() {
    try {
        final URL url = new URL("https://registry.hub.docker.com/v2/repositories/nanozoo/pangolin-v4/tags/");
        final URLConnection conn = url.openConnection();
        conn.connect();
        conn.getInputStream().close();
        return true;
    } catch (MalformedURLException e) {
        return false;
    } catch (IOException e) {
        return false;
    }
}

def internetcheck = DockernetIsAvailable()

if (params.update) {
println "\033[0;33mWarning: Running --update might not be poreCov compatible!\033[0m"
    if ( internetcheck.toString() == "true" ) { 
        tagname = 'https://registry.hub.docker.com/v2/repositories/nanozoo/pangolin-v4/tags/'.toURL().text.split(',"name":"')[1].split('","')[0]
        params.pangolindocker = "nanozoo/pangolin-v4:" + tagname
        println "\033[0;32mFound latest pangolin container, using: " + params.pangolindocker + " \033[0m" 

        tagname = 'https://registry.hub.docker.com/v2/repositories/nanozoo/nextclade3/tags/'.toURL().text.split(',"name":"')[1].split('","')[0]
        params.nextcladedocker = "nanozoo/nextclade3:" + tagname 
        println "\033[0;32mFound latest nextclade3 container, using: " + params.nextcladedocker + " \033[0m"
    } 
    if ( internetcheck.toString() == "false" ) { 
        println "\033[0;33mCould not find the latest pangolin container, trying: " + params.defaultpangolin + "\033[0m"
        params.pangolindocker = params.defaultpangolin 

        println "\033[0;33mCould not find the latest nextclade3 container, trying: " + params.defaultnextclade + "\033[0m"
        params.nextcladedocker = params.defaultnextclade 
    } 
}
else { params.pangolindocker = params.defaultpangolin ; params.nextcladedocker = params.defaultnextclade  }

if ( params.screen_reads && params.lcs_ucsc_update ){
    if ( internetcheck.toString() == "true" ) { 
        latest_version = 'https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt'.toURL().text.split('\\(')[1].split('\\)')[0]
        log.info "\033[0;32mFound latest UCSC version, using: " + latest_version + " \033[0m" 
        params.lcs_ucsc = latest_version
    }
    if ( internetcheck.toString() == "false" ) { 
        log.info "\033[0;33mCould not find the latest UCSC version, trying: " + params.lcs_ucsc_version + "\033[0m"
        params.lcs_ucsc = params.lcs_ucsc_version
    }
} else { params.lcs_ucsc = params.lcs_ucsc_version}


/************************** 
* Log-infos
**************************/

defaultMSG()
if ( params.fast5 || workflow.profile.contains('test_fast5') ) { basecalling() }
if (!params.fasta && !workflow.profile.contains('test_fasta')) { read_length() }
rki()

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

/************************** 
* MAIN WORKFLOW
**************************/

workflow {
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
                    reporterrorfast5 = basecalling_wf.out[0].join(samples_input_ch).ifEmpty{ exit 2, "Could not match barcode numbers from $params.samples to the read files, some typo?"} 
                    }
                else if (!params.samples) { fastq_from5_ch = basecalling_wf.out[0] }
            
            filtered_reads_ch = filter_fastq_by_length(fastq_from5_ch)
            noreadsatall = filtered_reads_ch.ifEmpty{ log.info "\033[0;33mNot enough reads in all samples, please investigate $params.output/$params.readqcdir\033[0m" }
            read_classification_wf(filtered_reads_ch)

            artic_ncov_wf(filtered_reads_ch, params.artic_normalize) 
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
                reporterrorfastq = fastq_input_raw_ch.join(samples_input_ch).ifEmpty{ exit 2, "Could not match barcode numbers from $params.samples to the read files, some typo?"} 
                }
                else if (!params.samples) { fastq_input_ch = fastq_input_raw_ch }

            read_qc_wf(fastq_input_ch)
            filtered_reads_ch = filter_fastq_by_length(fastq_input_ch)
            noreadsatall = filtered_reads_ch.ifEmpty{ log.info "\033[0;33mNot enough reads in all samples, please investigate $params.output/$params.readqcdir\033[0m" }
            read_classification_wf(filtered_reads_ch)

            // genome reconstruction with artic
            artic_ncov_wf(filtered_reads_ch, params.artic_normalize)
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

        determine_lineage_wf(fasta_input_ch)
        determine_mutations_wf(fasta_input_ch)
        genome_quality_wf(fasta_input_ch, reference_for_qc_input_ch)

    // 3. Specialised outputs (rki, json)
        rki_report_wf(genome_quality_wf.out[0], genome_quality_wf.out[1], extended_input_ch)

        if (params.samples) {
            create_json_entries_wf(determine_lineage_wf.out, genome_quality_wf.out[0], determine_mutations_wf.out)
        }

    // 4. Summary output
        if (params.fasta || workflow.profile.contains('test_fasta')) {
            taxonomic_read_classification_ch = Channel.from( ['deactivated', 'deactivated', 'deactivated'] ).collect()
            alignments_ch = Channel.from( ['deactivated'] )

        } else {
            taxonomic_read_classification_ch = read_classification_wf.out.kraken
            if (params.screen_reads) {
                if (params.lcs) {
                    read_screening_lsc_wf(filtered_reads_ch)
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
        create_summary_report_wf(determine_lineage_wf.out, genome_quality_wf.out[0], determine_mutations_wf.out,
                                taxonomic_read_classification_ch, alt_allele_ratio_ch, alignments_ch, samples_file_ch)

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
    .    
\033[0;33mUsage examples:${c_reset}

    nextflow run replikation/poreCov --update --fastq '*.fasta.gz' -r 1.3.0 -profile local,singularity

${c_yellow}Inputs (choose one):${c_reset}
    --fast5         one fast5 dir of a nanopore run containing multiple samples (barcoded);
                    to skip demultiplexing (no barcodes) add the flag [--single]
                    ${c_dim}[Basecalling + Genome reconstruction + Lineage + Reports]${c_reset}

    --fastq         one fastq or fastq.gz file per sample or
                    multiple file-samples: --fastq 'sample_*.fastq.gz'
                    ${c_dim}[Genome reconstruction + Lineage + Reports]${c_reset}

    --fastq_pass    the fastq_pass dir from the (guppy) bascalling
                    --fastq_pass 'fastq_pass/'
                    to skip demultiplexing (no barcodes) add the flag [--single]
                    ${c_dim}[Genome reconstruction + Lineage + Reports]${c_reset}

    --fasta         direct input of genomes - supports multi-fasta file(s) - can be gzip compressed (.gz)
                    ${c_dim}[Lineage + Reports]${c_reset}

${c_yellow}Workflow control (optional)${c_reset}
    --artic_normalize        Normalise down to moderate coverage to save runtime [default: $params.artic_normalize]
                                   ${c_dim}(after mapping and before variant calling in the ARTIC bioinformatics pipeline)
                                   Use `--artic_normalize False` to turn off this normalisation.${c_reset}
    --update                 Always try to use latest pangolin & nextclade release [default: $params.update]
    --samples                .csv input (header: Status,_id), renames barcodes (Status) by name (_id), e.g.:
                             Status,_id
                             barcode01,sample2011XY
                             BC02,thirdsample_run
    --list                   --fastq takes a csv file containing (new) sample names and paths to the fastq-files (no header).
                             Paths need to start with '/' or poreCov searches the files in the current working dir, e.g.:
                             sample1,/path_to_first_sample.fastq.gz
                             sample2,/path_to_second_sample.fastq.gz
    --extended               poreCov utilizes from --samples these additional headers:
                             Submitting_Lab,Isolation_Date,Seq_Reason,Sample_Type
    --screen_reads           Determines the Pangolineage of each individual read (takes time, needs --freyja and/or --lcs)
    --scorpio  Skip Scorpio in pangolin run [default: $params.scorpio]
                                  ${c_dim}From pangolin version 4, Scorpio overwrites Usher results which leads to many unassigned samples
                                  Can be turned on with --scorpio${c_reset}

${c_yellow}Parameters - Lineage detection on reads (see screen_reads, optional)${c_reset}
    --freyja    activate read-screening via freyja
    --freyja_update update freyja's barcode-db prior to running
    --lcs       activate read-screening via lcs
    --lcs_ucsc_version       Create marker table based on a specific UCSC SARS-CoV-2 tree (e.g. '2022-05-01'). Use 'predefined' 
                             to use the marker table from the repo (most probably not up-to-date) [default: $params.lcs_ucsc_version]
                                 ${c_dim}See https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2 for available trees.${c_reset}
    --lcs_ucsc_update        Use latest UCSC SARS-CoV-2 tree for marker table update. Overwrites --lcs_ucsc_version [default: $params.lcs_ucsc_update]
                                 ${c_dim}Automatically checks https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt${c_reset}
    --lcs_ucsc_downsampling  Downsample sequences when updating marker table to save resources. Use 'None' to turn off [default: $params.lcs_ucsc_downsampling]
                                 ${c_dim}Attention! Updating without downsampling needs a lot of resources in terms of memory and might fail.
                                 Consider downsampling or increase the memory for this process.${c_reset}
    --lcs_variant_groups     Provide path to custom variant groups table (TSV) for marker table update (requires --lcs_ucsc_update). Use 'default' 
                                 for predefined groups from repo (https://github.com/rki-mf1/LCS/blob/master/data/variant_groups.tsv) [default: $params.lcs_variant_groups]
    --lcs_cutoff             Plot linages above this threshold [default: $params.lcs_cutoff]     

${c_yellow}Parameters - Basecalling  (optional)${c_reset}
    --localguppy    use a native installation of guppy instead of a gpu-docker or gpu_singularity
    --guppy_cpu     use cpus instead of gpus for basecalling
    --one_end       removes the recommended "--require_barcodes_both_ends" from guppy demultiplexing
                    try this if to many barcodes are unclassified (beware - results might not be trustworthy)
    --guppy_model   guppy basecalling model [default: ${params.guppy_model}]
                    e.g. "dna_r9.4.1_450bps_hac.cfg" or "dna_r9.4.1_450bps_sup.cfg"

${c_yellow}Parameters - SARS-CoV-2 genome reconstruction (optional)${c_reset}
    --primerV               Supported primer variants or primer bed files - choose one [default: ${params.primerV}]
                                ${c_dim}ARTIC:${c_reset} V1.0.0, V2.0.0, V3.0.0, V4.0.0, V4.1.0, V.5.0.0, V.5.1.0, V.5.3.2 (length 400)
                                ${c_dim}NEB:${c_reset} VarSkipV1a, VarSkipV2, VarSkipV2b
                                ${c_dim}Other:${c_reset} V1200, V5.2.0 (lenght 1200) ${c_dim}(also known as midnight)${c_reset}
                                ${c_dim}Primer bed file:${c_reset} e.g. primers.bed  ${c_dim}See Readme for more help${c_reset}
    --schemeLength          primer scheme length, e.g. 400, 1200 [default: ${params.schemeLength}]
    --rapid                 rapid-barcoding-kit was used [default: ${params.rapid}]
    --minLength             min length filter raw reads [default: 100]
    --maxLength             max length filter raw reads [default: 700 (primer-scheme: V1-4, rapid); 1500 (primer-scheme: V1200, V5.2.0_1200)]
    --min_depth             nucleotides below min depth will be masked to "N" [default ${params.min_depth}]
    --clair3_model_dir      directory to look for clair3 model files [default: ${params.clair3_model_dir}]
    --clair3_model_name     clair3 model for the artic workflow [default: ${params.clair3_model_name}]

${c_yellow}Parameters - Genome quality control  (optional)${c_reset}
    --reference_for_qc      reference FASTA for consensus qc (optional, wuhan is provided by default)
    --seq_threshold         global pairwise ACGT sequence identity threshold [default: ${params.seq_threshold}]
    --n_threshold           consensus sequence N threshold [default: ${params.n_threshold}] 

${c_yellow}Options  (optional)${c_reset}
    --cores         amount of cores for a process (local use) [default: $params.cores]
    --max_cores     max amount of cores for poreCov to use (local use) [default: $params.max_cores]
    --memory        available memory [default: $params.memory]
    --output        name of the result folder [default: $params.output]
    --cachedir      defines the path where singularity images are cached
                    [default: $params.cachedir]
    --krakendb      provide a .tar.gz kraken database [default: auto downloads one]

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

       Note: The singularity profile automatically passes the following environment variables to the container. 
       to ensure execution on HPCs: HTTPS_PROXY, HTTP_PROXY, http_proxy, https_proxy, FTP_PROXY, ftp_proxy
    """.stripIndent()
}

def header(){
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    log.info """
________________________________________________________________________________
    
${c_green}poreCov${c_reset} | A Nextflow SARS-CoV-2 workflow for nanopore data
    """
}

def defaultMSG(){
    log.info """
    .
    \u001B[32mProfile:             $workflow.profile\033[0m
    \033[2mCurrent User:        $workflow.userName
    Nextflow-version:    $nextflow.version
    poreCov-version:     $workflow.revision
    \u001B[0m
    Pathing:
    \033[2mWorkdir location [-work-Dir]:
        $workflow.workDir
    Output dir [--output]: 
        $params.output
    Databases location [--databases]:
        $params.databases
    Singularity cache dir [--cachedir]: 
        $params.cachedir
    \u001B[1;30m______________________________________\033[0m
    Parameters:
    \033[2mClair3 model:         $params.clair3_model_name [--clair3_model_name]
    Clair3 model dir:     $params.clair3_model_dir
    Min depth nucleotide: $params.min_depth [--min_depth]
    Latest Pangolin/Nextclade?: $params.update [--update]
    CPUs to use:          $params.cores [--cores]
    Memory in GB:         $params.memory [--memory]\u001B[0m
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}

def basecalling() {
    log.info """
    Basecalling options:
    \033[2mUse local guppy?        $params.localguppy [--localguppy]  
    One end demultiplexing? $params.one_end [--one_end]
    Basecalling via CPUs?   $params.guppy_cpu [--guppy_cpu]
    Basecalling modell:     $params.guppy_model [--guppy_model]
    Rapid-barcode-kit?:     $params.rapid [--rapid]\u001B[0m
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}

def rki() {
    log.info """
    RKI output for german DESH upload:
    \033[2mOutput stored at:      $params.output/$params.rkidir  
    Min Identity to NC_045512.2: $params.seq_threshold [--seq_threshold]
    Min Depth used:        $params.min_depth [--min_depth]
       Min Depth should be 20 or more for RKI upload
    Proportion cutoff N:   $params.n_threshold [--n_threshold]\u001B[0m
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}

def read_length() {
    log_msg_read_min_length = params.minLength
    log_msg_read_max_length = params.maxLength

    if ( params.primerV.matches('V1200') || params.primerV.matches('V5.2.0_1200') || params.schemeLength == 1200) {

        if ( params.primerV.matches('V1200') || params.primerV.matches('V5.2.0_1200')){
            println "\033[0;33mWarning: Definition of primer scheme length via --primerV is deprecated, please use --schemeLength instead. Setting length to 1200.. \033[0m"
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
    \033[2mMin read-length set to: $log_msg_read_min_length [--minLength]
    Max read-length set to: $log_msg_read_max_length [--maxLength]\u001B[0m
    \u001B[1;30m______________________________________\033[0m
    """.stripIndent()
}
