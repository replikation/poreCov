process pangolin {
    label 'pangolin'
    publishDir "${params.output}/fasta/${name}/", mode: 'copy'
  input:
    tuple val(name), path(fasta)
  output:
    tuple val(name), path("lineage_report_${name}.csv") optional true
    path("lineage_report_${name}.csv") optional true
  script:
    """
    pangolin -t ${task.cpus} ${fasta}

    mv lineage_report.csv lineage_report_${name}.csv

    find . -name "*.csv" -size  0 -print -delete
    """
  }