process split_reference {
    label 'ubuntu'
    publishDir "${name}/fasta/", mode: 'copy'
  input:
    tuple val(name), path(fasta)
  output:
    path("${name}_contigs/*.fa")
  script:
    """
    split.sh ${name} ${fasta}
    """
  }
