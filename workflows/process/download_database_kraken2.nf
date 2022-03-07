process download_database_kraken2 {
    label "ubuntu"
        publishDir "${params.databases}/kraken2/", mode: 'copy'
        errorStrategy 'retry'
        maxRetries 1
    output:
        path("kraken.tar.gz")
    script:
    if (task.attempt.toInteger() == 1)
        """
        echo ${task.attempt}
        wget --no-check-certificate https://zenodo.org/record/6333909/files/GRCh38.p13_SC2_2022-03-01.tar.gz -O kraken.tar.gz
        """
    else if (task.attempt.toInteger() > 1)
        """
        echo ${task.attempt}
        wget --no-check-certificate https://osf.io/eprfq/download -O kraken.tar.gz
        """
    stub:
        """
        touch kraken.tar.gz
        """
}
