process snp_sites {
    label "snp_sites"
    input:
        tuple val(name), file(sites)
    output:
        tuple val(name), file("clean.core.aln")
    script:
        """
        snp-sites -c ${sites} > clean.core.aln
        """
}