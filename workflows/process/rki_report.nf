process rki_report {
    label "ubuntu"
    publishDir "${params.output}/${params.rkidir}/", mode: 'copy'
    input:
        path(pangolin_data)
        path(readme)
    output:
        tuple path("rki_report.csv"), path("${readme}")
    script:
        if (params.rki ==~ "[0-9]+")
        """
        rki_report_parser.sh ${params.rki}
        """
        else
        """
        rki_report_parser.sh "00000"
        """
}