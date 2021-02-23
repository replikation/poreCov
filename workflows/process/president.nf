process president {
    label "president"
    publishDir "${params.output}/${params.genomedir}/${name}", mode: 'copy',
        saveAs: { filename -> if (filename.endsWith("${name}_report.tsv")) "${name}_seq_ident_check.tsv" }
    input:
        tuple val(name), path(fasta), path(reference_fasta)
    output:
        tuple val(name), path("${name}_report.tsv"), path("${name}_valid.fasta"), emit: valid
        tuple val(name), path("${name}_report.tsv"), path("${name}_invalid.fasta"), emit: invalid
    script:
        """
        president -r ${reference_fasta} -t ${task.cpus} -q ${fasta} -x ${params.seq_threshold} -n ${params.n_threshold} -p . -f ${name}_
        """
}

