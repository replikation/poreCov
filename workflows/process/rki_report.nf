process rki_report {
    label "ubuntu"
    publishDir "${params.output}/${params.rkidir}/valid", mode: 'copy', pattern: "rki_valid_report.csv"
    publishDir "${params.output}/${params.rkidir}", mode: 'copy', pattern: "${readme}"
    input:
        path(president_data)
        path(readme)
    output:
        path("rki_valid_report.csv"), emit: report
        path("${readme}"), emit: readme
    script:
        """
        cat ${president_data} | sed '/(ignoring Ns)/d' | awk -F'\\t' '{print \$1 "\\t" \$8}' | sed '/\\False\\b/d' >> all.csv

        if [ -s all.csv ]; then
            rki_report_parser.sh all.csv rki_valid_report.csv
        else
            touch rki_valid_report.csv
        fi
        """
    stub:
        """
        touch rki_valid_report.csv ${readme}
        """
}

process rki_report_extended {
    label "ubuntu"
    publishDir "${params.output}/${params.rkidir}/valid", mode: 'copy', pattern: "rki_valid_report.csv"
    publishDir "${params.output}/${params.rkidir}", mode: 'copy', pattern: "${readme}"
    input:
        path(president_data)
        path(readme)
        path(extended_table)
    output:
        path("rki_valid_report.csv"), emit: report
        path("${readme}"), emit: readme
    script:
        """
        cat ${president_data} | sed '/(ignoring Ns)/d' | awk -F'\\t' '{print \$1 "\\t" \$8}' | sed '/\\False\\b/d' >> all.csv

        if [ -s all.csv ]; then
            rki_report_parser_extended.sh all.csv rki_valid_report.csv ${extended_table}
        else
            touch rki_valid_report.csv
        fi
        """
    stub:
        """
        touch rki_valid_report.csv ${readme}
        """
}

/* 
Adjust the awk if president changes:
        # awk field 1 = fastaname
        # awk field 8 = all qc true
*/