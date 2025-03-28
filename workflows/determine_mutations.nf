include { nextclade } from './process/nextclade'

workflow determine_mutations_wf {
    take: 
        fasta
        nextcladedocker
    main:
        nextclade(fasta, nextcladedocker)

    emit:
        nextclade.out
} 
