include { rki_report } from './process/rki_report'

workflow rki_report_wf {
    take: 
        pangolin_results  
    main:
        readme_pdf_ch = Channel
        .fromPath(workflow.projectDir + "/data/rki_report/Readme.pdf")

        rki_report(pangolin_results.map{it -> it[1]}.collect(), readme_pdf_ch)
        
    emit: 
        rki_report.out
} 