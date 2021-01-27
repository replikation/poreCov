include { rki_report } from './process/rki_report'

workflow rki_report_wf {
    take: 
        president_valid
        president_invalid
    main:
        readme_pdf_ch = Channel
        .fromPath(workflow.projectDir + "/data/rki_report/Readme.pdf")

        rki_report(president_valid.map{it -> it[1]}.collect(), readme_pdf_ch)

        // store valid genomes
        channel_tmp1 = president_valid.map{it -> it[2]}
            .splitText(by:100000000)
            .collectFile(name: 'valid_genomes.fasta', storeDir: params.output + "/" + params.rkidir +"/valid/")

        // store invalid genomes
        channel_tmp2 = president_invalid.map{it -> it[2]}
            .splitText(by:100000000)
            .collectFile(name: 'invalid_genomes.fasta', storeDir: params.output + "/" + params.rkidir +"/invalid/")

    emit: 
        rki_report.out.report
} 