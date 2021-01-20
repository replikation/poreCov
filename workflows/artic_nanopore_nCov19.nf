include { artic; artic_V1200 } from './process/artic.nf' 
include { bwa_samtools } from './process/bwa_samtools.nf'
include { coverage_plot } from './process/coverage_plot.nf'
include { filter_fastq_by_length } from './process/filter_fastq_by_length.nf'

workflow artic_ncov_wf {
    take:   
        fastq
    main: 

        // assembly
        if ( params.primerV.matches('V1200') ) {
            external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            artic_V1200(filter_fastq_by_length(fastq).combine(external_primer_schemes))

            assembly = artic_V1200.out.fasta
        }
        else {
            artic(filter_fastq_by_length(fastq))

            assembly = artic.out.fasta
        }

        // validate fasta
        coverage_plot(
            bwa_samtools(
                assembly.join(filter_fastq_by_length.out))[0])

    emit:   
        assembly
}