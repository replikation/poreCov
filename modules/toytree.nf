process toytree {

    publishDir "${params.output}/${clonal_complex}/tree/", mode: 'copy', pattern: "tree.svg"
    publishDir "${params.output}/${clonal_complex}/tree/", mode: 'copy', pattern: "tree.pdf"
    label 'toytree'
  input:
    tuple val(clonal_complex), val(refname), path(tree), path(metadata)
  output:
	  tuple path("tree.svg"), path("tree.pdf")
  script:
    """
    vis_tree_features.py --tree ${tree} --reroot ${refname} --data ${metadata}  --outfmt svg,pdf --prefix tree
    """
}



/*


*/
