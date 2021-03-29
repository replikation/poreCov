include { artic_medaka ; artic_nanopolish } from './process/artic.nf' 
include { bwa_samtools } from './process/bwa_samtools.nf'
include { coverage_plot } from './process/coverage_plot.nf'
include { filter_fastq_by_length } from './process/filter_fastq_by_length.nf'

workflow artic_ncov_wf {
    take:   
        fastq
    main: 

        // assembly
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            
            filter_fastq_by_length(fastq)

            artic_medaka(filter_fastq_by_length.out.combine(external_primer_schemes))

            assembly = artic_medaka.out.fasta

        // validate fasta
        coverage_plot(
            bwa_samtools(
                assembly.join(filter_fastq_by_length.out))[0])

        // error logging
        noreadsatall = filter_fastq_by_length.out.ifEmpty("\033[0;33mNot enough reads in all samples, please investigate $params.output/$params.readqcdir\033[0m").view()

    emit:   
        assembly
        filter_fastq_by_length.out
}

workflow artic_ncov_np_wf {
    take:   
        fastq
        fast5
        sequence_summaries
    main: 

        // assembly
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            artic_nanopolish(
                filter_fastq_by_length(fastq)
                    .combine(external_primer_schemes)
                    .combine(fast5.map{it -> it[1]})
                    .combine(sequence_summaries)
                    .map{it -> tuple(it[0],it[1],it[2],it[3],it[5])}
            )

            assembly = artic_nanopolish.out.fasta

        // validate fasta
        coverage_plot(
            bwa_samtools(
                assembly.join(filter_fastq_by_length.out))[0])

        // error logging
        noreadsatall = filter_fastq_by_length.out.ifEmpty("\033[0;33mNot enough reads in all samples, please investigate $params.output/$params.readqcdir\033[0m").view()

    emit:   
        assembly
        filter_fastq_by_length.out
}