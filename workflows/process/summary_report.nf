process summary_report {
        publishDir "${params.output}/", mode: 'copy'
        // re-use pangolin container for pandas dependency
        label 'pangolin'
    input:
        path(pangolin_results)
        path(president_results)
        path(nextclade_results)
    output:
	    path("poreCov_summary_report.html")
    script:
    """
    summary_report.py \
        -p ${pangolin_results} \
        -q ${president_results} \
        -n ${nextclade_results}
    """
}