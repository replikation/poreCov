process quality_genome_filter {
    label 'ubuntu'
    publishDir "${params.output}/fasta/${name}/", mode: 'copy', pattern: "error_report_*.txt"
  input:
    tuple val(name), path(fasta)
  output:
    tuple val(name), path("qc_${fasta}") optional true
    tuple val(name), path("error_report_*.txt") optional true
  script:
    """
    genome_integrity.sh ${fasta} ${params.maskBegin} ${params.maskEnd} ${params.rm_N_genome}

    if [ ! -f error_report_*.txt ]
      then
          mv ${fasta} qc_${fasta}         
    fi
    """
  }