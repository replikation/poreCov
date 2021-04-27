process pangolin {
    label 'pangolin'
    container = params.pangolindocker
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "lineage_report_${name}.csv"
  input:
    tuple val(name), path(fasta)
  output:
    tuple val(name), path("lineage_report_${name}.csv") optional true
  script:
    """
    pangolin -t ${task.cpus} ${fasta}

    mv lineage_report.csv lineage_report_${name}.csv
    
    find . -name "*.csv" -size  0 -print -delete
    """
  }