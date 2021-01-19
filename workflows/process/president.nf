process president {
    label "president"
    publishDir "${params.output}/fasta/${name}", mode: 'copy', pattern: '*.tsv'
    input:
        tuple val(name), path(fasta), path(reference_fasta)
    output:
        tuple val(name), path("${name}_seq_ident_check.tsv"), optional: true
    script:
        """
        president -r ${reference_fasta} -q ${fasta} -x ${params.threshold} -o ${name}_seq_ident_check.tsv
        """
}

