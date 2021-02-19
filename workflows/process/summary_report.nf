process summary_report {
        publishDir "${params.output}/", mode: 'copy'
        // re-use pangolin container for pandas dependency
        label 'pangolin'
    input:
        path(pangolin_results)
        path(president_results)
        path(nextclade_results)
        file(kraken2_results)
        path(version_config)
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
        echo '#sample,num_sarscov2,num_human' > kraken2_results.csv
        for KF in !{kraken2_results}; do
        echo -n "${KF%.kreport}," >> kraken2_results.csv
        NSARS=$(awk -v ORS= '$5=="2697049" {print $3 "," }' $KF)
        echo ${NSARS:-0} >> kraken2_results.csv
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
            -k kraken2_results.csv
        '''
}