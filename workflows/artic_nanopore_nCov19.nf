include { artic_medaka ; artic_nanopolish; artic_medaka_custom_bed; artic_nanopolish_custom_bed } from './process/artic.nf' 
include { covarplot; covarplot_custom_bed } from './process/covarplot.nf'

workflow artic_ncov_wf {
    take:   
        fastq
        normalise_threshold
    main: 

        // assembly with a primer bed file
        if (params.primerV.toString().contains(".bed")) {
            primerBed = Channel.fromPath(params.primerV, checkIfExists: true )
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )

            artic_medaka_custom_bed(fastq.combine(external_primer_schemes).combine(primerBed), normalise_threshold)
            assembly = artic_medaka_custom_bed.out.fasta
            binary_alignment = artic_medaka_custom_bed.out.fullbam
            trimmed_bam = artic_medaka_custom_bed.out.reference_bam
            vcf = artic_medaka_custom_bed.out.vcf
            failed_vcf = artic_medaka_custom_bed.out.vcf_fail
            primer_dir = artic_medaka_custom_bed.out.primer_dir

            // plot amplicon coverage
            covarplot_custom_bed(artic_medaka_custom_bed.out.covarplot.combine(primerBed))
        }
            
        // assembly via pre installed Primers
        else {
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            
            artic_medaka(fastq.combine(external_primer_schemes), normalise_threshold)
            assembly = artic_medaka.out.fasta
            binary_alignment = artic_medaka.out.fullbam
            trimmed_bam = artic_medaka.out.reference_bam
            vcf = artic_medaka.out.vcf
            failed_vcf = artic_medaka.out.vcf_fail
            primer_dir = Channel.empty()

            // plot amplicon coverage
            covarplot(artic_medaka.out.covarplot.combine(external_primer_schemes))
        }


        // error logging
        no_genomes_at_all = assembly.ifEmpty{ log.info "\033[0;33mCould not generate any genomes, please check your reads $params.output/$params.readqcdir\033[0m" }
        binary_alignment.ifEmpty{ log.info "\033[0;33mCould not generate any genomes, please check your reads $params.output/$params.readqcdir\033[0m" }

    emit:   
        assembly
        binary_alignment
        trimmed_bam
        vcf
        primer_dir
        failed_vcf
}

workflow artic_ncov_np_wf {
    take:   
        fastq
        fast5
        sequence_summaries
        normalise_threshold
    main: 

        // assembly
        if (params.primerV.toString().contains(".bed")) {
            primerBed = Channel.fromPath(params.primerV, checkIfExists: true )
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )

            artic_nanopolish_custom_bed(
                fastq
                    .combine(external_primer_schemes)
                    .combine(fast5.map{it -> it[1]})
                    .combine(sequence_summaries)
                    .combine(primerBed)
                    .map{it -> tuple(it[0],it[1],it[2],it[3],it[5],it[6])},
                normalise_threshold
        )

            assembly = artic_nanopolish_custom_bed.out.fasta
            binary_alignment = artic_nanopolish_custom_bed.out.fullbam

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
                        .map{it -> tuple(it[0],it[1],it[2],it[3],it[5])},
                    normalise_threshold
            )

            assembly = artic_nanopolish.out.fasta
            binary_alignment = artic_nanopolish.out.fullbam

            // plot amplicon coverage
            covarplot(artic_nanopolish.out.covarplot.combine(external_primer_schemes))
        }
        
    emit:   
        assembly
        binary_alignment
}
