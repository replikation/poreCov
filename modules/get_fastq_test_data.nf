process get_nanopore_fastq {
  	label 'template'
	output:
	path("SARSCoV2.fastq.gz") 
	script:
	"""
    wget --no-check-certificate https://osf.io/kf54a/download -O SARSCoV2.fastq.gz
	"""
}
