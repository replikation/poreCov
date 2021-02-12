include { summary_report } from './process/summary_report' 

workflow create_summary_report_wf {
    take: 
        pangolin
        president
        nextclade
    main:
        pangolin_results = pangolin.map {it -> it[1]}.collectFile(name: 'pangolin_results.csv', skip: 1, keepHeader: true)
        president_results = president.map {it -> it[1]}.collectFile(name: 'president_results.tsv', skip: 1, keepHeader: true)
        nextclade_results = nextclade.map {it -> it[1]}.collectFile(name: 'nextclade_results.tsv', skip: 1, keepHeader: true)

        summary_report(pangolin_results, president_results, nextclade_results)

} 
