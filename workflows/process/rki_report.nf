process rki_report {
    label "ubuntu"
    publishDir "${params.output}/${params.rkidir}/valid", mode: 'copy', pattern: "rki_valid_report.csv"
    publishDir "${params.output}/${params.rkidir}", mode: 'copy', pattern: "${readme}"
    input:
        path(pangolin_data)
        path(readme)
    output:
        path("rki_valid_report.csv"), emit: report
        path("${readme}"), emit: readme
    script:
        if (params.rki ==~ "[0-9]+")
        """
        rki_report_parser.sh ${params.rki} rki_valid_report.csv
        """
        else
        """
        rki_report_parser.sh "00000" rki_valid_report.csv
        """
}