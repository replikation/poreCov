process filter_fastq_by_length {
    label 'ubuntu'
  input:
    tuple val(name), path(reads) 
  output:
	  tuple val(name), path("${name}_filtered.fastq.gz") 
  script:
    """
    case "${reads}" in
      *.fastq.gz ) 
        zcat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= 400' | sed 's/\\t/\\n/g' |\
          awk -F"\\t" 'length(\$2)  <= 700' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
        ;;
      *.fastq)
        cat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= 400' | sed 's/\\t/\\n/g' |\
          awk -F"\\t" 'length(\$2)  <= 700' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
        ;;
    esac   
    """
}

/* Comments:
This is a super fast process to remove short reads.

it can take .fastq or .fastq.gz
*/