process guppy_gpu {
        maxForks 1
        if (!params.localguppy && workflow.profile.contains('docker') ) {
            container = 'nanozoo/guppy_gpu:4.2.2-1--5fc71df'
            containerOptions '--gpus all'
        }
        publishDir "${params.output}/fastq/", mode: 'copy'
    input:
        tuple val(name), path(dir)
    output:
        tuple val(name), path("*.fastq.gz"), emit: reads
        tuple val(name), path("fastq_tmp/*.txt"), emit: summary
    script:
        if (params.single)
        """
        guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i ${dir} -s fastq -x auto -r --trim_strategy dna -q 0

        find -L fastq -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
        
        mkdir -p fastq_tmp/
        cp fastq/*.txt fastq_tmp
        """
        else if (!params.single && !params.one_end)
        """
        guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i ${dir} -s fastq_tmp -x auto -r
        guppy_barcoder --require_barcodes_both_ends -i fastq_tmp -s fastq --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"

        for barcodes in fastq/barcode??; do
            find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
        done

        cp fastq/*.txt fastq_tmp
        """
        else if (!params.single && params.one_end)
        """
        guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -i ${dir} -s fastq_tmp -x auto -r
        guppy_barcoder -i fastq_tmp -s fastq --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"

        for barcodes in fastq/barcode??; do
            find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
        done

        cp fastq/*.txt fastq_tmp
        """
}
