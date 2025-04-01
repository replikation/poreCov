include { president } from './process/president' 
include { seqrs } from './process/seqrs' 

workflow genome_quality_wf {
    take:
        fasta
        reference
    main:
        president(fasta.combine(reference))

        primerbedfile = Channel.fromPath(workflow.projectDir + "/data/external_primer_schemes/nCoV-2019/", checkIfExists: true, type: 'dir' )
        seqrs(fasta.combine(primerbedfile))

    emit:   
        president_valid = president.out.valid
        president_invalid = president.out.invalid
}