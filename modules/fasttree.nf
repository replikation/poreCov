process fasttree {
    label "fasttree"
    publishDir "${name}/tree/", mode: 'copy', pattern: "${name}_clean.core.tree.nwk"
    input:
        tuple val(name), file(clean_core_alignment)
    output:
        tuple val(name), file(clean_core_alignment), file("${name}_clean.core.tree.nwk")
    script:
        """
        FastTree -gtr -nt ${clean_core_alignment} > ${name}_clean.core.tree.nwk
        """
}