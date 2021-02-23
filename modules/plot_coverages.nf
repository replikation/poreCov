
process plot_coverages {
  	label 'fastcov'
    input:
        path(alignment_files)
		path(index_files)
	output:
	    path('coverages.png')
	script:
	"""
    fastcov.py -l -o coverages.png ${alignment_files}
	"""
}
