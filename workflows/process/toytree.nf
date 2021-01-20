process toytree {
        publishDir "${params.output}/tree/", mode: 'copy'
        label 'toytree'
        //errorStrategy 'retry'
        maxRetries 1
    input:
        path(tree)
    output:
	    tuple path("tree*.svg"), path("tree*.pdf")
    script:
    if (task.attempt.toString() == '1')
    """
    render_tree.py --tree ${tree} --highlight ${params.highlight} --format 1
    render_tree_circle.py --tree ${tree} --highlight ${params.highlight} --format 1
    """
    else if (task.attempt.toString() == '2')
    """
    render_tree.py --tree ${tree} --highlight ${params.highlight} --format 0
    render_tree_circle.py --tree ${tree} --highlight ${params.highlight} --format 0
    """
}