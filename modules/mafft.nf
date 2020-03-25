process mafft {
    label "mafft"
    input:
        tuple val(name), path(sample)
        tuple val(reference_name), path(references)
    output:
        tuple val(name), path(sample), path("clean.full.aln")
    script:
        """
        cat ${sample} ${references} > combined.fasta
        mafft --thread ${task.cpus}--auto combined.fasta > clean.full.aln
        """
}

