process mafft {
    label "mafft"
    input:
        path(samples)
        tuple val(reference_name), path(references)
    output:
        path("clean.full.aln")
    script:
        """
        # merge sequence and add counter for uniq names
        cat ${samples} ${references} | tr -d "\\r" | tr "|" "_" | tr "/" "_" | tr " " "_" | awk '/^>/{\$0=\$0"_"(++i)}1' > combined.fasta
        mafft --thread ${task.cpus} combined.fasta > clean.full.aln
        """
}

