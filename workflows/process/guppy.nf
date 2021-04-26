process guppy_gpu {
    label 'guppy_gpu'
        if (!params.localguppy) {
            if (workflow.profile.contains('docker')) {
                container = 'nanozoo/guppy_gpu:4.4.1-1--a3fcea3'
                containerOptions '--gpus all'
            }
            else if (workflow.profile.contains('singularity')) {
                container = 'nanozoo/guppy_gpu:4.4.1-1--a3fcea3'
                containerOptions '--nv'
            }
            else if (workflow.profile.contains('ukj_cloud') || workflow.profile.contains('nanozoo')) {
            accelerator 2, type: 'nvidia-tesla-p100'
            container = 'nanozoo/guppy_gpu:4.4.1-1--a3fcea3'
            containerOptions '--gpus all'
            }
        }
        else if (params.localguppy) {
            if (workflow.profile.contains('ukj_cloud') || workflow.profile.contains('nanozoo')) {
                executor = "local"
                container = 'nanozoo/guppy_gpu:4.4.1-1--a3fcea3'
                containerOptions '--gpus all'
            }
            else {
                executor = "local"
            }
        }

        errorStrategy { if ( task.exitStatus == 127) { 'retry' ; exit 1, "Could not find the guppy basecaller"  }
                    else if (task.exitStatus == 255) { 'retry' ; exit 1, "nvidia docker toolkit not installed (correctly)?" }
                    else if (task.exitStatus == 125) { 'retry' ; exit 1, "nvidia cuda driver not found" } }

        publishDir "${params.output}/${params.readsdir}/", mode: 'copy'
    input:
        tuple val(name), path(dir)
    output:
        tuple val(name), path("*.fastq.gz"), emit: reads
        tuple val(name), path("fastq_tmp/*.txt"), emit: summary
    script:       
        if ( params.kit.contains('RBK') ) {
            guppy_arrangement_files = 'barcode_arrs_rbk4.cfg'
            barcoding_option = '  '
            }
        else {
            barcoding_option = '--require_barcodes_both_ends'
            guppy_arrangement_files = 'barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg'
            }
        if ( params.one_end ) {
            barcoding_option = '  '
            }
        if (params.single)
        """
        guppy_basecaller -c ${params.guppy_model} -i ${dir} -s fastq -x auto -r --trim_strategy dna -q 0

        find -L fastq -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
        
        mkdir -p fastq_tmp/
        cp fastq/*.txt fastq_tmp
        """
        else
        """
        guppy_basecaller -c ${params.guppy_model} -i ${dir} -s fastq_tmp -x auto -r
        guppy_barcoder -t ${task.cpus} ${barcoding_option} -i fastq_tmp -s fastq --arrangements_files "${guppy_arrangement_files}"

        for barcodes in fastq/barcode??; do
            find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
        done

        cp fastq/*.txt fastq_tmp
        """
}

process guppy_cpu {
        label 'guppy_cpu'
        if (!params.localguppy && workflow.profile.contains('docker') || workflow.profile.contains('singularity') ) {
            container = 'nanozoo/guppy_cpu:4.2.2-1--416f83d'
        }
        publishDir "${params.output}/${params.readsdir}/", mode: 'copy'
    input:
        tuple val(name), path(dir)
    output:
        tuple val(name), path("*.fastq.gz"), emit: reads
        tuple val(name), path("fastq_tmp/*.txt"), emit: summary
    script:
        if ( params.kit.contains('RBK') ) {
            guppy_arrangement_files = 'barcode_arrs_rbk4.cfg'
            barcoding_option = '  '
            }
        else {
            barcoding_option = '--require_barcodes_both_ends'
            guppy_arrangement_files = 'barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg barcode_arrs_nb96.cfg'
            }
        if ( params.one_end ) {
            barcoding_option = '  '
            }
        if (params.single)
        """
        guppy_basecaller -c ${params.guppy_model} -i ${dir} -s fastq  --num_callers ${task.cpus} --cpu_threads_per_caller 1 -r --trim_strategy dna -q 0

        find -L fastq -name '*.fastq' -exec cat {} +  | gzip > ${name}.fastq.gz
        
        mkdir -p fastq_tmp/
        cp fastq/*.txt fastq_tmp
        """
        else
        """
        guppy_basecaller -c ${params.guppy_model} -i ${dir} -s fastq_tmp  --num_callers ${task.cpus} --cpu_threads_per_caller 1 -r
        guppy_barcoder -t ${task.cpus} ${barcoding_option} -i fastq_tmp -s fastq --arrangements_files "${guppy_arrangement_files}"

        for barcodes in fastq/barcode??; do
            find -L \${barcodes} -name '*.fastq' -exec cat {} + | gzip > \${barcodes##*/}.fastq.gz
        done

        cp fastq/*.txt fastq_tmp
        """
}
