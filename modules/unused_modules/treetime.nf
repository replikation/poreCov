process treetime {
    label "treetime"
    publishDir "${params.output}/${clonal_complex}/", mode: 'copy', pattern: "${clonal_complex}_tree_corrected.newick" 
    input:
        tuple val(name), file(clean_core_alignment), file(tree), file(dates)
    output:
        tuple val(clonal_complex), val(ref_name), file("${clonal_complex}_tree_corrected.newick")
    script:
        """
        treetime --keep-root --max-iter 5 --gtr infer --tip-labels \
            --aln ${clean_core_alignment} --tree ${tree} --dates ${dates} --outdir .

        # convert nexus to newick (from bin/)
        nexus2newick.py timetree.nexus ${clonal_complex}_tree.nwk

        #sed 's/'Reference'/'${ref_name}'/' ${clonal_complex}_tree.nwk  > ${clonal_complex}_tree_corrected.newick
        """
}