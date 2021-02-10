include { json_report } from './process/json_report' 

workflow create_json_entries_wf {
    take: 
        pangolin
        president
    main:
        merged_ch = pangolin.join(president.map { it -> tuple(it[0], it[1])})

        json_report(merged_ch)

        if (params.fasta) {
            // sample csv merge
        }
        

} 
