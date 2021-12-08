process add_aainsertions {
    label 'fastcov'
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}_aaIns_clade.tsv"
    input:
        tuple val(name), path(nextclade_result)
    output:
        tuple val(name), path("${name}_aaIns_clade.tsv")
    script:
    """
    convert_insertions_nt2aa.py --nextclade_results ${nextclade_result} --output ${name}_aaIns_clade.tsv
    """
    stub:
    """
    touch ${name}_aaIns_clade.tsv
    """
}
