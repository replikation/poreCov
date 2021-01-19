
process get_fasta {
  	container = 'nanozoo/template:3.8--d089809'
	storeDir "tmp_input_test_files/nanopore_fastq" 
	output:
	path("SARSCoV2.fasta") 
	script:
	"""
    wget https://osf.io/87bc9/download -O SARSCoV2.fasta
	"""
}
