include { summary_report } from './process/summary_report' 

workflow create_summary_report_wf {
    take: 
        pangolin
        president
        nextclade
        kraken2
    main:
        pangolin_results = pangolin.map {it -> it[1]}.collectFile(name: 'pangolin_results.csv', skip: 1, keepHeader: true)
        president_results = president.map {it -> it[1]}.collectFile(name: 'president_results.tsv', skip: 1, keepHeader: true)
        nextclade_results = nextclade.map {it -> it[1]}.collectFile(name: 'nextclade_results.tsv', skip: 1, keepHeader: true)
        kraken2_results = kraken2.map {it -> it[2]}.collect()
        version_ch = Channel.fromPath(workflow.projectDir + "/configs/container.config")

        summary_report(pangolin_results, president_results, nextclade_results, kraken2_results, version_ch)

} 
