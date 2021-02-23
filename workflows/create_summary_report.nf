include { summary_report } from './process/summary_report'
include { plot_coverages } from '../modules/plot_coverages.nf'

workflow create_summary_report_wf {
    take: 
        pangolin
        president
        nextclade
        kraken2
        alignments
    main:
        version_ch = Channel.fromPath(workflow.projectDir + "/configs/container.config")

        pangolin_results = pangolin.map {it -> it[1]}.collectFile(name: 'pangolin_results.csv', skip: 1, keepHeader: true)
        president_results = president.map {it -> it[1]}.collectFile(name: 'president_results.tsv', skip: 1, keepHeader: true)
        nextclade_results = nextclade.map {it -> it[1]}.collectFile(name: 'nextclade_results.tsv', skip: 1, keepHeader: true)
        kraken2_results = kraken2.map {it -> it[2]}.collect()
        alignment_files = alignments.map {it -> it[0]}.collect()
        if (params.fasta || workflow.profile.contains('test_fasta')) {
            coverage_plot = Channel.from( ['deactivated'] )
        } else {
            coverage_plot = plot_coverages(alignments.map{it -> it[0]}.collect(), alignments.map{it -> it[1]}.collect())
        }

        summary_report(version_ch, pangolin_results, president_results, nextclade_results, kraken2_results, coverage_plot)

} 
