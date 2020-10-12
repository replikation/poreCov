process artic {
        label 'artic'
        publishDir "${params.output}/fasta/${name}/", mode: 'copy'
    input:
        tuple val(name), path(reads)
    output:
        tuple val(name), path("*.consensus.fasta")
    script:
        """
        artic minion --medaka --normalise 200 --threads ${task.cpus} --scheme-directory /artic-ncov2019/primer_schemes \
            --read-file ${reads} nCoV-2019/${params.primerV} ${name}
        """
}

process artic_V1200 {
        label 'artic'
        publishDir "${params.output}/fasta/${name}/", mode: 'copy'
    input:
        tuple val(name), path(reads), path(external_scheme)
    output:
        tuple val(name), path("*.consensus.fasta")
    script:
        if ( params.primerV.matches('V1200') )
        """
        artic minion --medaka --normalise 200 --threads ${task.cpus} --scheme-directory ${external_scheme} \
            --read-file ${reads} nCoV-2019/${params.primerV} ${name}
        """

}