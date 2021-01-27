include { president } from './process/president' 

workflow genome_quality_wf {
    take:
        fasta
        reference
    main:
        president(fasta.combine(reference))
    emit:   
        president.out.valid
        president.out.invalid
}