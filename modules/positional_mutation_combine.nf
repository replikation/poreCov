process positional_mutation_combiner {
    label 'fastcov'
    publishDir "${params.output}/${params.lineagedir}/", mode: 'copy', pattern: "positional_mutation_summary.csv"
    input:
        path(filePaths)
    output:
        path("positional_mutation_summary.csv")
    script:
        """
        positional_mutation_combine.py -f ${filePaths}
        """
        stub:
        """
        touch positional_mutation_summary.csv
        """
}