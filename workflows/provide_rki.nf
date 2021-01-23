include { rki_report } from './process/rki_report'

workflow rki_report_wf {
    take: 
        pangolin_results
        president
    main:
        readme_pdf_ch = Channel
        .fromPath(workflow.projectDir + "/data/rki_report/Readme.pdf")

        rki_report(pangolin_results.map{it -> it[1]}.collect(), readme_pdf_ch)
    
        //#### waiting here on president #37 to be fixed
        valid_reports = president.map { it -> it[1] }
            .splitCsv(header: true, sep: '\t')
            .filter{ it -> it[1] == 'True' } // only using entries with a "true" value
            //.map { it -> it [0]} 
            // modify rki rki_report.out
            // write to file
        
        invalid_reports = president.map { it -> it[1] }
            .splitCsv(header: true, sep: '\t')
            .filter{ it -> it[1] == 'False' } // only using entries with a "false" value
            //.map { it -> it [0]} 
            // write to file rki_report.out

        //####

        // store valid genomes
        president.map{it -> it[2]}
                .flatten()
                .filter{ it =~ '_valid.fasta' }
                .collectFile(name: 'valid_genomes.fasta', storeDir: params.output + "/" + params.rkidir +"/valid/")

        // store invalid genomes
        president.map{it -> it[2]}
            .flatten()
            .filter{ it =~ '_invalid.fasta' }
            .collectFile(name: 'invalid_genomes.fasta', storeDir: params.output + "/" + params.rkidir +"/invalid/")


    emit: 
        rki_report.out.report
} 