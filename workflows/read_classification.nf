include { kraken2 } from './process/kraken2.nf' 
include { krona } from './process/krona.nf' 
include { download_database_kraken2 } from './process/download_database_kraken2.nf'

workflow read_classification_wf {
    take:   
        fastq
    main: 
        download_database_kraken2()

        // trimming primer away is missing here (samclip macht das mit softclipped bases, aber hard coded nicht nein. seqtk?)
        kraken2(fastq, download_database_kraken2.out)

        // visuals
        krona(kraken2.out)

    emit:   
        kraken2.out
}