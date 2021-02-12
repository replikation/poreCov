process nextclade {
    label 'nextclade'
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}_clade.tsv"
    input:
        tuple val(name), path(consensus)
    output:
        tuple val(name), path("${name}_clade.tsv")
    """
    nextclade --input-fasta ${consensus} --output-tsv tmp.tsv
    cat tmp.tsv | tr -d "\r" > ${name}_clade.tsv
    """
}