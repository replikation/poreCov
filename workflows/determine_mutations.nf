include { nextclade } from './process/nextclade'
include { add_aainsertions } from '../modules/add_aainsertions.nf'

workflow determine_mutations_wf {
    take: 
        fasta  
    main:
        nextclade(fasta)
        add_aainsertions(nextclade.out)

    emit:
        add_aainsertions.out
} 
