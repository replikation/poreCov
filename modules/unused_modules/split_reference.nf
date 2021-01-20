process split_reference {
    label 'ubuntu'
  input:
    tuple val(name), path(fasta)
  output:
    path("${name}_contigs/*.fa")
  script:
    """
    split.sh ${name} ${fasta}
    """
  }
