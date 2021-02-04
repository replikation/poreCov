include { kraken2 } from './process/kraken2.nf' 
include { download_database_kraken2 } from './process/download_database_kraken2.nf'

workflow read_classification_wf {
    take:   
        fastq
    main: 
        download_database_kraken2()

        // trimming primer away is missing here
        kraken2(fastq, download_database_kraken2.out)

        // visuals are missing
        // also there is a kreport output similar to centrifuge to salvage some visuals
        // parse to krona (https://github.com/marbl/Krona/issues/117)
        // parse to sankey (https://github.com/EBI-Metagenomics/emg-viral-pipeline/blob/master/nextflow/modules/sankey.nf)
        // kraken needs also a database input flag to circumvent some issues

    emit:   
        kraken2.out
}