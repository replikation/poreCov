process president {
    label "president"
    publishDir "${params.output}/${params.genomedir}/${name}", mode: 'copy',
        saveAs: { filename -> if (filename.endsWith("${name}_report.tsv")) "${name}_seq_ident_check.tsv" }
    input:
        tuple val(name), path(fasta), path(reference_fasta)
    output:
        tuple val(name), path("output/${name}_report.tsv"), path("output/*_valid.fasta"), emit: valid
        tuple val(name), path("output/${name}_report.tsv"), path("output/*_invalid.fasta"), emit: invalid
    script:
        """
        president -r ${reference_fasta} -t ${task.cpus} -q ${fasta} -x ${params.threshold} -p output/${name}

        ## workarounds for current president error on invalid fasta ##
        if [ ! -f output/${name}_report.tsv ]; then
            echo "ID	Valid	ACGT Nucleotide identity	ACGT Nucleotide identity (ignoring Ns)	ACGT Nucleotide identity (ignoring non-ACGTNs)	Ambiguous Bases	Query Length	Query #ACGT	Query #IUPAC-ACGT	Query #non-IUPAC:	aligned	passed_initial_qc	Date	reference_length	reference	query" > output/${name}_report.tsv
            echo "False" >> output/${name}_report.tsv
            cat ${fasta} > output/${name}_invalid.fasta
            touch output/${name}_valid.fasta
        fi
        ## workaround end ##
        """
}

