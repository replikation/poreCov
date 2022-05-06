include { artic_medaka ; artic_nanopolish; artic_medaka_custom_bed; artic_nanopolish_custom_bed; artic_scheme_validation } from './process/artic.nf' 
include { covarplot; covarplot_custom_bed } from './process/covarplot.nf'

workflow artic_ncov_wf {
    take:   
        fastq
    main: 

        // assembly with a primer bed file
        if (params.primerV.toString().contains(".bed")) {
            primerBed = Channel.fromPath(params.primerV, checkIfExists: true )
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )

            artic_scheme_validation(primerBed)
            artic_scheme_validation.out.view()
            //if (artic_scheme_validation.out.head.contains('error')) {
            //    exit 10 "Primer.bed-file contains an error. Please check the report-file in '${params.output}/X.Pipeline-runinfo' for further information"
            //}

            artic_medaka_custom_bed(fastq.combine(external_primer_schemes).combine(primerBed))
            assembly = artic_medaka_custom_bed.out.fasta

            // plot amplicon coverage
            covarplot_custom_bed(artic_medaka_custom_bed.out.covarplot.combine(primerBed))
        }
            
        // assembly via pre installed Primers
        else {
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            
            artic_medaka(fastq.combine(external_primer_schemes))
            assembly = artic_medaka.out.fasta

            // plot amplicon coverage
            covarplot(artic_medaka.out.covarplot.combine(external_primer_schemes))
        }


        // error logging
        no_genomes_at_all = assembly.ifEmpty{ log.info "\033[0;33mCould not generate any genomes, please check your reads $params.output/$params.readqcdir\033[0m" }

    emit:   
        assembly
}

workflow artic_ncov_np_wf {
    take:   
        fastq
        fast5
        sequence_summaries
    main: 

        // assembly
        if (params.primerV.toString().contains(".bed")) {
            primerBed = Channel.fromPath(params.primerV, checkIfExists: true )
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            
            artic_scheme_validation(primerBed)
            //if (artic_scheme_validation.out.head.contains("error")) {
            //    exit 10 "Primer.bed-file contains an error. Please check the report-file in ${params.output}/X.Pipeline-runinfo for further information"
            //}

            artic_nanopolish_custom_bed(
                fastq
                    .combine(external_primer_schemes)
                    .combine(fast5.map{it -> it[1]})
                    .combine(sequence_summaries)
                    .combine(primerBed)
                    .map{it -> tuple(it[0],it[1],it[2],it[3],it[5],it[6])}
        )

            assembly = artic_nanopolish_custom_bed.out.fasta

            // plot amplicon coverage
            covarplot_custom_bed(artic_nanopolish_custom_bed.out.covarplot.combine(primerBed))
        }


        // assembly via pre installed Primers
        else {
        external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            artic_nanopolish(
                    fastq
                        .combine(external_primer_schemes)
                        .combine(fast5.map{it -> it[1]})
                        .combine(sequence_summaries)
                        .map{it -> tuple(it[0],it[1],it[2],it[3],it[5])}
            )

            assembly = artic_nanopolish.out.fasta

            // plot amplicon coverage
            covarplot(artic_nanopolish.out.covarplot.combine(external_primer_schemes))
        }
        
    emit:   
        assembly
}
