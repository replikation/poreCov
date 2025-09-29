
process plot_coverages {
  	label 'fastcov'

    input:
        path(alignment_files)
		path(index_files)
	output:
	    path("coverages_*.png")
	script:
	"""
	output_name=coverages_\$(echo ${alignment_files} | tr ' ' '_').png
	output_name_no_ext=\${output_name//.bam/}
	output_name=\$output_name_no_ext

	if [ \${#output_name} -gt 255 ]; then
		echo "Output file name exceeds filename character limit"
		output_name="\${output_name:0:240}.png"
	fi

    fastcov.py -l -o \$output_name ${alignment_files}
	fastcov.py -l -p NC_045512.2:21563-25385 -o \$output_name ${alignment_files}
	"""
	stub:
	"""
	touch coverages_1.png
	"""
}
