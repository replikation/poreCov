include { guppy_gpu } from './process/guppy'
include { guppy_cpu } from './process/guppy'
include { pycoqc } from './process/pycoqc'

workflow basecalling_wf {
    take: 
        dir_input  
    main:
        
        if (params.guppy_cpu) { 
            guppy_cpu(dir_input) 
            guppy_basecalls = guppy_cpu.out.reads
            guppy_summary = guppy_cpu.out.summary
        }
        else if (!params.guppy_cpu) { 
            guppy_gpu(dir_input)
            guppy_basecalls = guppy_gpu.out.reads
            guppy_summary = guppy_gpu.out.summary
        }

        // nanopore sequencing summary
        pycoqc(guppy_summary)

        // adjust channels to a clean val(name), path(reads) channel
        if (params.single) { fastq_channel = guppy_basecalls }
        else { fastq_channel = guppy_basecalls
                            .map { it -> it[1] }
                            .flatten()
                            .map { it -> [ it.simpleName, it ] }
        }

    emit:
        fastq_channel
} 