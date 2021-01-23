process president {
    label "president"
    publishDir "${params.output}/${params.genomedir}/${name}", mode: 'copy',
        saveAs: { filename -> if (filename.endsWith(".tsv")) "${name}_seq_ident_check.tsv" }
    input:
        tuple val(name), path(fasta), path(reference_fasta)
    output:
        tuple val(name), path("output/${name}_report.tsv"), path("output/*valid.fasta"), optional: true
    script:
        """
        president -r ${reference_fasta} -t ${task.cpus} -q ${fasta} -x ${params.threshold} -p output/${name}
        """
}

