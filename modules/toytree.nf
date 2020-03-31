process toytree {
        publishDir "${params.output}/${name}/tree/", mode: 'copy'
        label 'toytree'
        errorStrategy 'retry'
        maxRetries 1
    input:
        tuple val(name), path(tree)
    output:
	    tuple path("tree.svg"), path("tree.pdf")
    script:
    if (task.attempt.toString() == '1')
    """
    cleanname=\$(printf "${name}" | cut -f 1 -d ".")

    render_tree.py --tree ${tree} --highlight \${cleanname} --format 1
    """
    else if (task.attempt.toString() == '2')
    """
    cleanname=\$(printf "${name}" | cut -f 1 -d ".")

    render_tree.py --tree ${tree} --highlight \${cleanname} --format 0
    """
}