include { nextclade } from './process/nextclade'

workflow determine_mutations_wf {
    take: 
        fasta  
    main:
        nextclade(fasta, "/data/sars-cov-2_MN908947")

    emit:
        nextclade.out
} 
