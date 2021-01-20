include { create_database } from './process/create_database'

workflow build_database_wf {
    main:
        fasta_DB = Channel.fromPath( workflow.projectDir + "/database/ena_*.fasta" , checkIfExists: true)
        text_DB = Channel.fromPath( workflow.projectDir + "/database/ena_*.txt", checkIfExists: true)
    
        create_database(fasta_DB, text_DB)
    emit:
        create_database.out[0]
        create_database.out[1]
}