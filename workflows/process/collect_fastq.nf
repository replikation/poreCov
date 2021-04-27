process collect_fastq {
        label 'demultiplex'
    input:
        tuple val(name), path(dir)
    output:
        tuple val(name), path("*.fastq.gz"), emit: reads
    script:
        if (params.single)
        """
        find -L ${dir} -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
        """
        else if (!params.single )
        """
        BARCODE_DIRS=\$(find -L ${dir} -name "barcode??" -type d)
        
        if [ -z "\${BARCODE_DIRS}" ]; then 
            guppy_barcoder -t ${task.cpus} -r --require_barcodes_both_ends -i ${dir} -s fastq_porecov --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg"

            for barcodes in fastq_porecov/barcode??; do
                find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
            done
        else
            for barcodes in \${BARCODE_DIRS}; do
                find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
                find -L \${barcodes} -name '*.fastq.gz' -exec zcat {} + | gzip >> \${barcodes##*/}.fastq.gz
            done
        fi
        """
}
