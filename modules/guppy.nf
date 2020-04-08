process guppy_gpu {
        maxForks 1
        container = 'nanozoo/guppy_gpu:3.4.4-1--3dd30a4'
        containerOptions '--gpus all'
        publishDir "${params.output}/fastq/", mode: 'copy'
    input:
        tuple val(name), path(dir)
    output:
        tuple val(name), path("*.fastq.gz")
    script:
        if (!params.barcodes)
        """
        guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i ${dir} -s fastq -x auto -r --trim_strategy dna -q 0

        find -L fastq -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
        """
        else if (params.barcodes)
        """
        guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i ${dir} -s ${name} -x auto -r
        guppy_barcoder  -i ${name} -s fastq --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"

        #--require_barcodes_both_ends

        for barcodes in fastq/barcode??; do
            find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
        done
        """
}
