process artic {
    label 'artic'  
  input:
    tuple val(name), path(reads)
  output:
    tuple val(name), path(read), path("${name}.fastq.consensus.fasta")
  script:
    """
    artic minion --medaka --normalise 200 --threads ${task.cpus} --scheme-directory /artic-ncov2019/primer_schemes \
    --read-file ${reads} nCoV-2019/${params.primerV} ${name}
    """
  }

  /*

      artic-ncov2019/primer-schemes/nCoV-2019/V1/nCoV-2019.reference.fasta
      artic-ncov2019/primer_schemes/nCoV-2019/V1


    #mkdir -p fastq_dir
    #mv ${reads} fastq_dir
    #artic gather --min-length 400 --max-length 700 --prefix run_name --directory fastq_dir
    #artic gather --min-length 400 --max-length 700 --prefix run_name --directory /path/to/reads

    #artic demultiplex --threads 4 run_name_pass.fastq

  */