include { collect_fastq } from './process/collect_fastq'

workflow collect_fastq_wf {
    take: 
        fastq_dir  
    main:
        collect_fastq(fastq_dir)

        if (params.single) { fastq_channel = collect_fastq.out }
        else { fastq_channel = collect_fastq.out
                            .map { it -> it[1] }
                            .flatten()
                            .map { it -> [ it.simpleName, it ] }
        }

    emit: fastq_channel
} 