process cat_fastq {
    label 'ubuntu'
  input:
    tuple val(name), path(fastq_dir) 
  output:
	  tuple val(name), path("${name}_filtered.fastq.gz") 
  script:
    """
    find ${fastq_dir} -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
    """
}

/* Comments:
This is a super fast process to remove short reads.

it can take .fastq or .fastq.gz
*/