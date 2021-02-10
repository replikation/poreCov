process download_database_kraken2 {
    label "ubuntu"
        publishDir "${params.databases}/kraken2/", mode: 'copy'
    output:
        path("kraken.tar.gz")
    script:
        """
        wget https://zenodo.org/record/3854856/files/GRCh38.p13_GBcovid19-2020-05-22.tar.gz?download=1 -O kraken.tar.gz
        """
}
