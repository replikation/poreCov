include { summary_report } from './process/summary_report' 

workflow create_summary_report_wf {
    take: 
        pangolin
        president
        nextclade
    main:
        merged_ch = pangolin.join(president.map { it -> tuple(it[0], it[1])}.join(nextclade))

        summary_report(merged_ch)

} 
