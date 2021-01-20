include { nanoplot } from './process/nanoplot'

workflow read_qc_wf {
    take: 
        fastq  
    main:
        nanoplot(fastq)
} 