process pangolin {
    label 'pangolin'
    container { params.pangolindocker }
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "lineage_report_${name}.csv"
  input:
    tuple val(name), path(fasta)
  output:
    tuple val(name), path("lineage_report_${name}.csv"), optional: true
  script:
    def args = params.scorpio ? '' : '--skip-scorpio'
    """
    pangolin ${args} -t ${task.cpus} ${fasta}

    mv lineage_report.csv lineage_report_${name}.csv
    
    find . -name "*.csv" -size  0 -print -delete
    """
    stub:
    """
    echo "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note" > lineage_report_${name}.csv
    echo "barcode13_ARTIC_medaka,B.1.177,,,,,,PANGO-v1.2.12,3.0.5,2021-06-05,v1.2.12,passed_qc,Assigned from designation hash." >> lineage_report_${name}.csv
    """
  }
