process toytree {
    publishDir "${params.output}/${name}/tree/", mode: 'copy'
    label 'toytree'
  input:
    tuple val(name), path(tree)
  output:
	  tuple path("tree.svg"), path("tree.pdf")
    
    //def clean_name = name.split('.')
  script:
    """
    cleanname=\$(printf "${name}" | cut -f 1 -d ".")

    render_tree.py --tree ${tree} --highlight \${cleanname}
    """
}