include { devider; freyja_covariants } from './process/devider.nf'
// adjust if freyja_covariants gets moved to freyja.nf

workflow ww_cryptic_lineages_wf {
    take:
        fastq
        reference
        bam_and_bai
        
    main:
        devider(fastq.combine(reference)) 
        freyja_covariants(bam_and_bai.combine(reference))

    emit:
        devider.out
        freyja_covariants.out

}