process nextclade {
    label 'nextclade'
    container { params.nextcladedocker }
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}_clade.tsv"
    input:
        tuple val(name), path(consensus)
    output:
        tuple val(name), path("${name}_clade.tsv")
    script:
    """
    nextclade run -j ${task.cpus} --input-dataset /data/sars-cov-2_MN908947 --output-tsv tmp.tsv ${consensus}

    cat tmp.tsv | tr -d "\r" > ${name}_clade.tsv
    """
    stub:
    """
    touch ${name}_clade.tsv
    """
}