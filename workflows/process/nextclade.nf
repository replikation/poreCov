process nextclade {
    label 'nextclade'
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}_clade.tsv"
    input:
        tuple val(name), path(consensus)
        val(reference_file)
    output:
        tuple val(name), path("${name}_clade.tsv")
    script:
    """
    nextclade run --input-fasta ${consensus} --input-dataset ${reference_file} --output-tsv tmp.tsv
    cat tmp.tsv | tr -d "\r" > ${name}_clade.tsv
    """
    stub:
    """
    touch ${name}_clade.tsv
    """
}