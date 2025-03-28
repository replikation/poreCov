process summary_report {

    publishDir "${params.output}/", mode: 'copy'
    label 'fastcov'
    
    input:
        path(version_config)
        path(variants_table)
        path(pangolin_results)
        path(president_results)
        path(nextclade_results)
        file(kraken2_results)
        file(mixed_sites_results)
        file(coverage_plots)
        file(samples_table)
        val(nextcladedocker)
    output:
	    path("poreCov_summary_report_*.html")
        path("poreCov_summary_report_*.xlsx")
        path("poreCov_summary_report_*.tsv")

    script:
        guppyused = ((params.fast5 && !(params.fastq || params.fastq_pass)) || workflow.profile.contains('test_fast5'))

        """
        echo 'sample,num_unclassified,num_sarscov2,num_human' > kraken2_results.csv
        for KF in ${kraken2_results}; do
        NUNCLASS=\$(awk -v ORS= '\$5=="0" {print \$3}' \$KF)
        NSARS=\$(awk -v ORS= '\$5=="2697049" {print \$3}' \$KF)
        NHUM=\$(awk '\$5=="9606" {print \$3}' \$KF)
        echo "\${KF%.kreport},\${NUNCLASS:-0},\${NSARS:-0},\${NHUM:-0}" >> kraken2_results.csv
        done

        summary_report.py \
            -v ${version_config} \
            --variants_table ${variants_table} \
            --porecov_version ${workflow.revision}:${workflow.commitId}:${workflow.scriptId} \
            --nextclade_docker ${nextcladedocker} \
            --guppy_used ${guppyused} \
            --guppy_model ${params.guppy_model} \
            --medaka_model ${params.medaka_model} \
            --nf_commandline '${workflow.commandLine}' \
            --primer ${params.primerV} \
            -p ${pangolin_results} \
            -q ${president_results} \
            -n ${nextclade_results} \
            -k kraken2_results.csv \
            -m ${mixed_sites_results} \
            -c \$(echo ${coverage_plots} | tr ' ' ',') \
            -s ${samples_table}
        """
    stub:
        """
        touch poreCov_summary_report_1.html poreCov_summary_report_1.xlsx poreCov_summary_report_1.tsv
        """
}

process summary_report_default {

    publishDir "${params.output}/", mode: 'copy'
    label 'fastcov'
    
    input:
        path(version_config)
        path(variants_table)
        path(pangolin_results)
        path(president_results)
        path(nextclade_results)
        path(kraken2_results)
        path(mixed_sites_results)
        path(coverage_plots)
        val(nextcladedocker)
    output:
	    path("poreCov_summary_report_*.html")
        path("poreCov_summary_report_*.xlsx")
        path("poreCov_summary_report_*.tsv")

    script:
        guppyused = ((params.fast5 && !(params.fastq || params.fastq_pass)) || workflow.profile.contains('test_fast5'))

        """
        echo 'sample,num_unclassified,num_sarscov2,num_human' > kraken2_results.csv
        for KF in ${kraken2_results}; do
        NUNCLASS=\$(awk -v ORS= '\$5=="0" {print \$3}' \$KF)
        NSARS=\$(awk -v ORS= '\$5=="2697049" {print \$3}' \$KF)
        NHUM=\$(awk '\$5=="9606" {print \$3}' \$KF)
        echo "\${KF%.kreport},\${NUNCLASS:-0},\${NSARS:-0},\${NHUM:-0}" >> kraken2_results.csv
        done

        summary_report.py \
            -v ${version_config} \
            --variants_table ${variants_table} \
            --porecov_version ${workflow.revision}:${workflow.commitId}:${workflow.scriptId} \
            --guppy_used ${guppyused} \
            --guppy_model ${params.guppy_model} \
            --medaka_model ${params.medaka_model} \
            --nf_commandline '${workflow.commandLine}' \
            --nextclade_docker ${nextcladedocker} \
            --primer ${params.primerV} \
            -p ${pangolin_results} \
            -q ${president_results} \
            -n ${nextclade_results} \
            -k kraken2_results.csv \
            -m ${mixed_sites_results} \
            -c \$(echo !{coverage_plots} | tr ' ' ',') 
        """
    stub:
        """
        touch poreCov_summary_report_1.html poreCov_summary_report_1.xlsx poreCov_summary_report_1.tsv
        """
}

process summary_report_fasta {
        publishDir "${params.output}/", mode: 'copy'
        label 'fastcov'
    input:
        path(version_config)
        path(variants_table)
        path(pangolin_results)
        path(president_results)
        path(nextclade_results)
        val(nextcladedocker)
    output:
	    path("poreCov_summary_report_*.html")
        path("poreCov_summary_report_*.xlsx")
        path("poreCov_summary_report_*.tsv")

    script:          
        """
        summary_report.py \
            -v ${version_config} \
            --variants_table ${variants_table} \
            --porecov_version ${workflow.revision}:${workflow.commitId}:${workflow.scriptId} \
            --nf_commandline '${workflow.commandLine}' \
            --nextclade_docker ${nextcladedocker} \
            -p ${pangolin_results} \
            -q ${president_results} \
            -n ${nextclade_results} \
            -s "deactivated"
        """
    stub:
        """
        touch poreCov_summary_report_1.html poreCov_summary_report_1.xlsx poreCov_summary_report_1.tsv
        """
}
