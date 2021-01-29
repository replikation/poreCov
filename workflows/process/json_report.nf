process json_report {
        publishDir "${params.output}/${params.jsondir}/", mode: 'copy'
        label 'ubuntu'
    input:
        tuple val(name), path(pangolin_result), path(president_result)
    output:
	    tuple val(name), path("*.json")
    script:
    """
    json_parser.sh -i ${name} -p ${pangolin_result} -q ${president_result} -r ${params.rki} -v ${params.primerV}
    """

}