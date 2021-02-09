process kraken2 {
        label 'kraken2'
        publishDir "${params.output}/${params.readqcdir}/1.read_classification", mode: 'copy'
    input:
        tuple val(name), path(reads)
        path(database)
  	output:
    	tuple val(name), path("${name}.kraken.out"), path("${name}.kreport")
  	script:
    """
    mkdir -p kraken_db && tar xzf ${database} -C kraken_db --strip-components 1

    # mask possible primer regions
    zcat ${reads} | sed '1b ; s/..............................\$/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/ ; n ' |\
        sed '1b ; s/^............................../NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/ ; n ' > masked_reads.fastq
    kraken2 --db kraken_db --threads ${task.cpus} --output ${name}.kraken.out --report ${name}.kreport masked_reads.fastq
    """
  }