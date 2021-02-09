include { nextclade } from './process/nextclade'

workflow determine_mutations_wf {
    take: 
        fasta  
    main:
        nextclade(fasta)

    emit:
        nextclade.out
} 
