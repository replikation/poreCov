
process get_fasta {
  	label 'template'
	output:
	path("SARSCoV2.fasta") 
	script:
	"""
    wget --no-check-certificate https://osf.io/87bc9/download -O SARSCoV2.fasta
	"""
}
