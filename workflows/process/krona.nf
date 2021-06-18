process krona {
        label 'krona'
        publishDir "${params.output}/${params.readqcdir}/1.read_classification", mode: 'copy'
    input:
        tuple val(name), path(kraken2), path(kreport)
  	output:
    	tuple val(name), file("${name}_krona.html")
  	script:
    """
    cat ${kreport} | cut -f 3,5 > file.krona
    ktImportTaxonomy file.krona -m 1
    mv *.html ${name}_krona.html
    """
    stub:
    """
    touch ${name}_krona.html
    """
}