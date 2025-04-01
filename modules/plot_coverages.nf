
process plot_coverages {
  	label 'fastcov'

    input:
        path(alignment_files)
		path(index_files)
	output:
	    path("coverages_*.png")
	script:
	"""
    fastcov.py -l -o coverages_\$(echo ${alignment_files} | tr ' ' '_').png ${alignment_files}
	fastcov.py -l -p NC_045512.2:21563-25385 -o coverages_spike_\$(echo ${alignment_files} | tr ' ' '_').png ${alignment_files}
	"""
	stub:
	"""
	touch coverages_1.png
	"""
}
