process json_report {
        publishDir "${params.output}/${params.jsondir}/", mode: 'copy'
        label 'fastcov'
    input:
        tuple val(name), path(pangolin_result), path(president_result), path(nextclade_result)
    output:
	    tuple val(name), path("*.json")
    script:
    if ((params.fasta || workflow.profile.contains('test_fasta')))
    """
    json_parser.py -i ${name} \
        -p ${pangolin_result} \
        -n ${nextclade_result} \
        -q ${president_result} \
    """
    else if (!params.fasta )
    """
    json_parser.py -i ${name} \
        -p ${pangolin_result} \
        -n ${nextclade_result} \
        -q ${president_result} \
        -v ${params.primerV}
    """

}