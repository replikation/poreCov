process json_report {
        publishDir "${params.output}/${params.jsondir}/", mode: 'copy'
        label 'ubuntu'
    input:
        tuple val(name), path(pangolin_result), path(president_result), path(nextclade_result)
    output:
	    tuple val(name), path("*.json")
    script:
    if ((params.fasta || workflow.profile.contains('test_fasta')) && params.rki)
    """
    json_parser.sh -i ${name} \
        -p ${pangolin_result} \
        -n ${nextclade_result} \
        -q ${president_result} \
        -r ${params.rki}
    """
    else if ((params.fasta || workflow.profile.contains('test_fasta')) && !params.rki)
    """
    json_parser.sh -i ${name} \
        -p ${pangolin_result} \
        -n ${nextclade_result} \
        -q ${president_result} 
    """
    else if (!params.fasta && params.rki)
    """
    json_parser.sh -i ${name} \
        -p ${pangolin_result} \
        -n ${nextclade_result} \
        -q ${president_result} \
        -v ${params.primerV} \
        -r ${params.rki}
    """
    else if (!params.fasta && !params.rki)
    """
    json_parser.sh -i ${name} \
        -p ${pangolin_result} \
        -n ${nextclade_result} \
        -q ${president_result} \
        -v ${params.primerV}
    """

}