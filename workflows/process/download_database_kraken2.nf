process download_database_kraken2 {
    label "ubuntu"
        publishDir "${params.databases}/kraken2/", mode: 'copy'
        errorStrategy 'retry'
        maxRetries 1
    output:
        path("kraken.tar.gz")
    script:
    if (task.attempt = 1)
        """
        wget --no-check-certificate https://zenodo.org/record/4534746/files/GRCh38.p13_SC2_2021-02-08.tar.gz -O kraken.tar.gz
        """
    if (task.attempt > 1)
        """
        wget --no-check-certificate https://osf.io/eprfq/download -O kraken.tar.gz
        """
    stub:
        """
        touch kraken.tar.gz
        """
}
