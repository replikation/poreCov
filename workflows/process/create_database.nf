process create_database {
    label "ubuntu"
    if (params.cloudProcess) { 
        publishDir "${params.databases}/references/", mode: 'copy' 
    }
    else { 
        storeDir "${params.databases}/references/" 
    } 
    input:
        path(fasta)
        path(embl)
    output:
        path("clean.references.fasta")
        path("clean_data.tsv")
    script:
        """
        create_references.sh ${fasta} ${embl}
        """
}