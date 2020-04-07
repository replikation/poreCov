process fasttree {
    label "fasttree"
    publishDir "${params.output}/tree/", mode: 'copy', pattern: "clean.core.tree.nwk"
    input:
        path(clean_core_alignment)
    output:
        tuple path(clean_core_alignment), path("clean.core.tree.nwk")
    script:
        """
        FastTree -gtr -nt ${clean_core_alignment} > clean.core.tree.nwk
        """
}