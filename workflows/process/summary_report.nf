process summary_report {
        publishDir "${params.output}/", mode: 'copy'
        label 'ubuntu'
    input:
        tuple val(name), path(pangolin_result), path(president_result), path(nextclade_result)
    output:
	    tuple val(name), path("*.html")
    script:
    """
    summary_report.py -i ${name} \
        -p ${pangolin_result} \
        -n ${nextclade_result} \
        -q ${president_result}
    """
}