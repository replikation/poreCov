process summary_report {
        publishDir "${params.output}/", mode: 'copy'
        label 'fastcov'
    input:
        path(version_config)
        path(pangolin_results)
        path(president_results)
        path(nextclade_results)
        file(kraken2_results)
        file(coverage_plots)
    output:
	    path("poreCov_summary_report_*.html")

    shell:
    if (params.fasta || workflow.profile.contains('test_fasta'))            
        '''
        summary_report.py \
            -v !{version_config} \
            --porecov_version !{workflow.revision}:!{workflow.commitId}:!{workflow.scriptId} \
            -p !{pangolin_results} \
            -q !{president_results} \
            -n !{nextclade_results}
        '''
    else
        '''
        echo '#sample,num_unclassified,num_sarscov2,num_human' > kraken2_results.csv
        for KF in !{kraken2_results}; do
        echo -n "${KF%.kreport}," >> kraken2_results.csv
        NUNCLASS=$(awk -v ORS= '$5=="0" {print $3 "," }' $KF)
        echo -n ${NUNCLASS:-0} >> kraken2_results.csv
        NSARS=$(awk -v ORS= '$5=="2697049" {print $3 "," }' $KF)
        echo -n ${NSARS:-0} >> kraken2_results.csv
        NHUM=$(awk '$5=="9606" {print $3}' $KF)
        echo ${NHUM:-0} >> kraken2_results.csv
        done

        summary_report.py \
            -v !{version_config} \
            --porecov_version !{workflow.revision}:!{workflow.commitId}:!{workflow.scriptId} \
            --primer !{params.primerV} \
            -p !{pangolin_results} \
            -q !{president_results} \
            -n !{nextclade_results} \
            -k kraken2_results.csv \
            -c $(echo !{coverage_plots} | tr ' ' ',')
        '''
}