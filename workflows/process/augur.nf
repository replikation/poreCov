/*
aligns sequences to a given reference, which is removed afterwards 
An alignmend guide so to say -> ${align_reference}
*/
process augur_align {
        label 'augur'
    input:
        path(sample)
        path(references)
		path(align_reference)
  	output:
    	path("aligned.fasta")
  	script:
    """
	# removing possible windows returns at line endings
	# cat ${sample} ${references} | tr -d "\\r" | tr "|" "_" | tr "/" "_" | tr " " "_" | awk '/^>/{\$0=\$0"_"(++i)}1' > combined.fasta
	cat ${sample} ${references} | tr -d "\\r" | tr ":|/ ,'" "_" > combined.fasta

    augur align \
        --sequences combined.fasta \
        --reference-sequence ${align_reference} \
        --output aligned.fasta \
        --fill-gaps \
        --remove-reference
    """
}


process augur_tree {
        label 'augur'
    input:
        path(alignment)
  	output:
    	tuple path(alignment), path("tree_raw.nwk")
  	script:
    """
	augur tree \
		--alignment ${alignment} \
		--output tree_raw.nwk
    """
}

process augur_tree_refine {
        label 'augur'
    input:
        tuple path(alignment), path(tree)
		path(metadata)
  	output:
    	path("tree_refined.nwk")
  	script:
    """
        augur refine \
            --tree ${tree} \
            --alignment ${alignment} \
            --metadata ${metadata} \
            --timetree \
            --output-tree tree_refined.nwk \
            --output-node-data nodes.json
    """
}

/*

//config/sars-like-cov_reference.gb

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
        """
    


rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.align.output.alignment
    output:
        alignment = "results/masked.fasta"
    params:
        mask_from_beginning = 15,
        mask_from_end = 15,
        mask_sites = 18460
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment}
        """



  */