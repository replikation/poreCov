include { summary_report; summary_report_fasta; summary_report_default } from './process/summary_report'
include { plot_coverages } from '../modules/plot_coverages.nf'
include { get_variants_classification } from '../modules/get_variants_classification.nf'
include { lcs_plot as lcs_plot; lcs_plot as lcs_plot_control  } from './process/lcs_sc2'

workflow create_summary_report_wf {
    take: 
        pangolin
        president
        nextclade
        kraken2
        lcs
        alignments
        samples_table

        main:
        version_ch = Channel.fromPath(workflow.projectDir + "/configs/container.config")
        variants_table_ch = get_variants_classification()

        pangolin_results = pangolin.map {it -> it[1]}.collectFile(name: 'pangolin_results.csv', skip: 1, keepHeader: true)
        president_results = president.map {it -> it[1]}.collectFile(name: 'president_results.tsv', skip: 1, keepHeader: true)
        nextclade_results = nextclade.map {it -> it[1]}.collectFile(name: 'nextclade_results.tsv', skip: 1, keepHeader: true)
        lcs_results = lcs.map {it -> it[1]}.collectFile(name: 'lcs_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.lineagedir}/")

        alignment_files = alignments.map {it -> it[0]}.collect()
        if (params.fasta || workflow.profile.contains('test_fasta')) {
            
            summary_report_fasta(version_ch, variants_table_ch, pangolin_results, president_results, nextclade_results)

        } else {
            kraken2_results = kraken2.map {it -> it[2]}.collect()
            // sort by sample name, group in lists of 6, collect the grouped plots
            coverage_plots = plot_coverages(alignments.map{it -> it[0]}.toSortedList({ a, b -> a.simpleName <=> b.simpleName }).flatten().collate(6), \
                                            alignments.map{it -> it[1]}.toSortedList({ a, b -> a.simpleName <=> b.simpleName }).flatten().collate(6)).collect()

            if (params.samples) { summary_report(version_ch, variants_table_ch, pangolin_results, president_results, nextclade_results, kraken2_results, coverage_plots, samples_table) }
            else { summary_report_default(version_ch, variants_table_ch, pangolin_results, president_results, nextclade_results, kraken2_results, coverage_plots) }
            
        }

        if (params.screen_reads){
            lcs_plot(lcs.filter({ !it[0].contains("negativ") }).map {it -> it[1]}.collectFile(name: 'lcs_results.tsv', skip: 1, keepHeader: true).map{ it -> ['samples', it] }, params.screen_reads_plot_cutoff)

            lcs_plot_control(lcs.filter({ it[0].contains("negativ") }), params.screen_reads_plot_cutoff)
        }

} 
