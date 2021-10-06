
process get_fasta {
  	container = 'nanozoo/template:3.8--d089809'
	output:
	path("SARSCoV2.fasta") 
	script:
	"""
    wget --no-check-certificate https://osf.io/87bc9/download -O SARSCoV2.fasta
	"""
}
