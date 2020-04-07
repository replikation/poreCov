process snp_sites {
    label "snp_sites"
    input:
        file(sites)
    output:
        file("clean.core.aln")
    script:
        """
        snp-sites -c ${sites} > clean.core.aln
        """
}