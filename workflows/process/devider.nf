process devider {
        label 'devider'
        publishDir "${params.output}/${params.lineagedir}/devider", mode: 'copy', pattern: "*_majority_vote_haplotypes.fasta" // should this be published or only intermediate for later?
        // probably change output dir

    input:
        tuple val(name), path(reads), path(reference)//, path(bam_file), path(vcf)
    output:
        tuple val(name), path("${name}_majority_vote_haplotypes.fasta"), emit: haplotypes
        // add other stuff if relevant

    script:
    // just run devider haplotyping
    /*"""
    devider -b ${bam_file} -v ${vcf} -r ${reference} -o ${name}_devider_output \
        -t ${task.cpus} --preset nanopore-r9
    """*/ // didn't work when I tried it separately (because it doesnt recognize chromosome name or smth)
    
    // run devider pipeline (does mapping, variant calling, haplotyping)
    """
    run_devider_pipeline -i ${reads} -r ${reference} -o ${name}_devider_output
    mv ${name}_devider_output/majority_vote_haplotypes.fasta ${name}_majority_vote_haplotypes.fasta
    """ 
    /* different presets - choose one or make dynamic??
        old-long-reads for ~90% identity long reads
        nanopore-r9 for ~95% (default)
        nanopore-r10 for ~98%
        hi-fi for true high fidelity reads (> 99.9% accuracy).
    */

}


process freyja_covariants {
        label 'freyja'
        publishDir "${params.output}/${params.lineagedir}/devider", mode: 'copy', pattern: "*"
        // probably change output dir
        // maybe move to freyja.nf for clarity?

    input:
        tuple val(name), path(bam_file), path(bam_index), path(reference)
        // needs .bai in same folder as .bam

    output:
        tuple val(name), path("${name}_freyja_covariants.tsv"), emit: covariants

    script:
    """
    freyja covariants --ref-genome ${reference} --output ${name}_freyja_covariants.tsv  ${bam_file} 
    """
    // could include .gff with --annot --> AA notes in result + nt mutation
    //  --threads ${task.cpus} gets error "no such option as --threads" but it exists according to --help and the wiki ??

}