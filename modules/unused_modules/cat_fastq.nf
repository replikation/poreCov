process cat_fastq {
        label 'ubuntu'
        echo true
    input:
        tuple val(name), path(fastq_dir) 
    output:
	    tuple val(name), path("*.fastq.gz") 
    script:
    if (!params.barcodes)
    """
    find -L ${fastq_dir} -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
    """
    else if (params.barcodes)
    """
    barcode_dirs=\$(find -L ${fastq_dir} -name 'barcode??')
    echo \$barcode_dirs
    for barcodes in \$barcode_dirs; do
        find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
    done
    """
}