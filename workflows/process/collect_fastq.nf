process collect_fastq {
        label 'demultiplex'
    input:
        tuple val(name), path(dir)
    output:
        tuple val(name), path("*.fastq.gz"), emit: reads
    script:
        if (params.rapid) {
            guppy_arrangement_files = 'barcode_arrs_rbk4.cfg barcode_arrs_rbk096.cfg'
            barcoding_option = '  '
            }
        else {
            guppy_arrangement_files = 'barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg'
            barcoding_option = '--require_barcodes_both_ends'
            }
        if (params.one_end) {
            barcoding_option = '  '
            }
        if (params.single)
        """
        find -L ${dir} -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
        find -L ${dir} -name '*.fastq.gz' -exec zcat {} + | gzip >> ${name}.fastq.gz
        """
        else if (!params.single)
        """
        BARCODE_DIRS=\$(find -L ${dir} -name "barcode??" -type d)
        
        if [ -z "\${BARCODE_DIRS}" ]; then 
            guppy_barcoder -t ${task.cpus} -r ${barcoding_option} -i ${dir} -s fastq_porecov --arrangements_files "${guppy_arrangement_files}"

            for barcodes in fastq_porecov/barcode??; do
                find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
            done
        else
            for barcodes in \${BARCODE_DIRS}; do
                find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip >> \${barcodes##*/}.fastq.gz
                find -L \${barcodes} -name '*.fastq.gz' -exec zcat {} + | gzip >> \${barcodes##*/}.fastq.gz
            done
        fi
        """
}
