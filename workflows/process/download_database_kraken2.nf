process download_database_kraken2 {
    label "ubuntu"
        publishDir "${params.databases}/kraken2/", mode: 'copy'
        errorStrategy 'retry'
        maxRetries 1
    output:
        path("kraken.tar.gz")
    script:
    if (task.attempt.toString() == '1')
        """
        wget https://zenodo.org/record/3854856/files/GRCh38.p13_GBcovid19-2020-05-22.tar.gz?download=1 -O kraken.tar.gz
        """
    if (task.attempt.toString() == '2')
        """
        wget https://osf.io/qxmeh/download -O kraken.tar.gz
        """
}
