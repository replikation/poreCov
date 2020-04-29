process mask_alignment {
        label 'augur'
    input:
      path(alignment)
  	output:
    	path("alignment_masked.fasta")
  	script:
    """
    mask_alignment.py \
        --alignment ${alignment} \
        --mask-from-beginning ${params.maskBegin} \
        --mask-from-end ${params.maskEnd} \
        --output alignment_masked.fasta
    """
  }



/*
--mask-sites  (sites to mask)
*/
