include { kraken2 } from './process/kraken2.nf' 
include { krona } from './process/krona.nf' 
include { download_database_kraken2 } from './process/download_database_kraken2.nf'
include { lcs_sc2 } from './process/lcs_sc2' 

workflow read_classification_wf {
    take:   
        fastq
    main: 

        // database download
        if (params.krakendb) { kraken_db = file("${params.krakendb}") }
        else  { download_database_kraken2(); kraken_db = download_database_kraken2.out } 

        // classification
        kraken2(fastq, kraken_db)

        // visuals
        krona(kraken2.out)

        // calculate mixed/ pooled samples using LCS, https://github.com/rvalieris/LCS
        if (params.screen_reads) {
            lcs_sc2(fastq)
            lcs_output = lcs_sc2.out
        } else {
            lcs_output = Channel.empty()
        }
        
    emit:   
        kraken = kraken2.out
        lcs = lcs_output
}
