include { artic_medaka ; artic_nanopolish } from './process/artic.nf' 
include { covarplot } from './process/covarplot.nf'

workflow artic_ncov_wf {
    take:   
        fastq
    main: 

        // assembly
        external_primer_schemes = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes", checkIfExists: true, type: 'dir' )
            
        artic_medaka(fastq.combine(external_primer_schemes))

        assembly = artic_medaka.out.fasta

        // plot amplicon coverage
        covarplot(
            artic_medaka.out.covarplot.combine(external_primer_schemes)
        )

        // error logging
        noreadsatall = filter_fastq_by_length.out.ifEmpty{ log.info "\033[0;33mNot enough reads in all samples, please investigate $params.output/$params.readqcdir\033[0m" }
        nogenomesatall = artic_medaka.out.fasta.ifEmpty{ log.info "\033[0;33mCould not generate any genomes, please check your reads $params.output/$params.readqcdir\033[0m" }

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
        covarplot(
            artic_nanopolish.out.covarplot.combine(external_primer_schemes)
        )
        
    emit:   
        assembly
}
