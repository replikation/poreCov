include { guppy_gpu } from './process/guppy'
include { pycoqc } from './process/pycoqc'

workflow basecalling_wf {
    take: 
        dir_input  
    main:
        
        guppy_gpu(dir_input)
        pycoqc(guppy_gpu.out.summary)

        if (params.single) { fastq_channel = guppy_gpu.out.reads }

        else { fastq_channel = guppy_gpu.out.reads
                            .map { it -> it[1] }
                            .flatten()
                            .map { it -> [ it.simpleName, it ] }
            }

    emit:
        fastq_channel
} 