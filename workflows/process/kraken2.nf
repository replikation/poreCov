process kraken2 {
        label 'kraken2'
        publishDir "${params.output}/${params.readqcdir}/read_classification", mode: 'copy'
    input:
        tuple val(name), path(reads)
        path(database)
  	output:
    	tuple val(name), path("${name}.kraken.out"), path("${name}.kreport")
  	script:
    """
    mkdir -p kraken_db && tar xzf ${database} -C kraken_db --strip-components 1
    kraken2 --db kraken_db --threads ${task.cpus} --gzip-compressed  --output ${name}.kraken.out --report ${name}.kreport ${reads}
    """
  }