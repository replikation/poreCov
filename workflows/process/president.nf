process president {
    label "president"
    publishDir "${params.output}/fasta/${name}", mode: 'copy', pattern: '*.tsv'
    input:
        tuple val(name), path(fasta), path(reference_fasta)
    output:
        tuple val(name), path("${name}_seq_ident_check.tsv"), optional: true
    script:
        """
        president -r ${reference_fasta} -t ${task.cpus} -q ${fasta} -x ${params.threshold} -p output/${name}
        cp output/${name}_report.tsv ${name}_seq_ident_check.tsv
        """
}

