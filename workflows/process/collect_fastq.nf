process collect_fastq {
        label 'ubuntu'
    input:
        tuple val(name), path(dir)
    output:
        tuple val(name), path("*.fastq.gz"), emit: reads 
    script:
        if (params.single)
        """
        find -L fastq -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
        """
        else if (!params.single )
        """
        for barcodes in ${dir}/barcode??; do
            find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
        done
        find . -name "*.fastq.gz" -type 'f' -size -1500k -delete
        """
}
