process positional_mutation_scan {
    label 'fastcov'
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}_positional_mutation_search.csv"
    input:
        tuple val(name), path(fasta), path(clade_tsv)
        path(mutation_list)
    output:
        tuple val(name), path("*positional_mutation_status.csv")
    script:
        """
        positional_mutation_scan.py -f ${fasta} -m ${mutation_list} -n ${name} -t ${clade_tsv}
        """
        stub:
        """
        touch ${name}_positional_mutation_status.csv
        """
}