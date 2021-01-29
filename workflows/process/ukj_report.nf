process ukj_report {
    publishDir "${params.output}/", mode: 'copy'

    input:
      tuple val(name), path(pangolin_result), path(president_result)
    output:
      path("*.json")
    script:
    """
    bash ukj_reporter.sh -i ${name} -p ${pangolin_result} -q ${president_result} -r ${params.rki} -v ${params.primerV}
    """
}