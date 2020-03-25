process gubbins {
    label "gubbins"
    input:
        tuple val(name), file(alignment)
    output:
        tuple val(name), file("gubbins.filtered_polymorphic_sites.fasta")
    script:
        """
        run_gubbins.py -p gubbins ${alignment}
        """
}

