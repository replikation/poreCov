process mafft {
    label "mafft"
    input:
        tuple val(name), path(sample)
        tuple val(reference_name), path(references)
    output:
        tuple val(name), path("clean.full.aln")
    script:
        """
        # merge sequence and add counter for uniq names
        cat ${sample} ${references} | tr -d "\\r" | tr "|" "_" | tr "/" "_" | tr " " "_" | awk '/^>/{\$0=\$0"_"(++i)}1' > combined.fasta
        mafft --thread ${task.cpus} combined.fasta > clean.full.aln
        """
}

