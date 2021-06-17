
process plot_coverages {
  	label 'fastcov'
    input:
        path(alignment_files)
		path(index_files)
	output:
	    path("coverages_*.png")
	shell:
	'''
    fastcov.py -l -o coverages_$(echo !{alignment_files} | tr ' ' '_').png !{alignment_files}
	'''
	stub:
	"""
	touch coverages_1.png
	"""
}
