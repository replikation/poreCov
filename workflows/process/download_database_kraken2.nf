process download_database_kraken2 {
    label "ubuntu"
    if (params.cloudProcess) { 
        publishDir "${params.databases}/kraken2/", mode: 'copy' 
    }
    else { 
        storeDir "${params.databases}/kraken2/" 
    } 
    output:
        path("GRCh38.p13_GBcovid19-2020-05-22.tar.gz")
    script:
        """
        wget https://zenodo.org/record/3854856/files/GRCh38.p13_GBcovid19-2020-05-22.tar.gz?download=1 -O GRCh38.p13_GBcovid19-2020-05-22.tar.gz
        """
}
