
process get_fast5 {
  	label 'template'
	output:
	path("fast5") 
	script:
	"""
    wget --no-check-certificate https://osf.io/vxd7f/download -O SARSCoV2.tar.gz
    tar xzf SARSCoV2.tar.gz 
	"""
}
