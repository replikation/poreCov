
process get_fast5 {
  	container = 'nanozoo/template:3.8--d089809'
	output:
	path("fast5") 
	script:
	"""
    wget https://osf.io/vxd7f/download -O SARSCoV2.tar.gz
    tar xzf SARSCoV2.tar.gz 
	"""
}
