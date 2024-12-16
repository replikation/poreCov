process hercules {
        label 'hercules'
        errorStrategy 'ignore'
        maxRetries 1
        publishDir "${params.output}/${params.lineagedir}/${name}/hercules", mode: 'copy', pattern: "*"

    input:
        tuple val(name), path(reads)
    output:
        tuple val(name), path("")

    script:
        """
            mkdir Data
            cp *.fastq* Data/
            /home/docker/CommonFiles/WWAnalysis.sh
        """
    stub:
        """
        """
}

