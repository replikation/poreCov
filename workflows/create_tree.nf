include { augur_align; augur_tree; augur_tree_refine } from './process/augur'
include { mask_alignment } from './process/mask_alignment'
include { quality_genome_filter } from './process/quality_genome_filter'


workflow create_tree_wf {
    take: 
        fasta       // the nCov fasta (own samples or reconstructed here)
        references  // multiple references to compare against
        metadata    // tsv file of meta data  strain country date
    main:

        align_reference = Channel.fromPath( workflow.projectDir + "/data/reference_nCov19/MN908947.gb", checkIfExists: true)

        quality_genome_filter(fasta)

        // from [val, file] to [files]
        collect_fasta = quality_genome_filter.out[0].map{ it -> it [1]}
                                                    .collect()

        augur_tree(
            mask_alignment(
                augur_align(collect_fasta, references, align_reference)))

        augur_tree_refine(augur_tree.out, metadata)

    emit:
        augur_tree_refine.out
}