
process get_variants_classification {
  	container = 'nanozoo/template:3.8--d089809'
	publishDir "${params.output}/", mode: 'copy'
	cache false
	// errorStrategy 'ignore'
	// this is set globally anyway
	output:
	path("SARSCoV2_variants_*.csv")
	shell:
	'''
	DATE=`date +"%Y-%m-%d--%H-%M-%S"`
    wget --no-check-certificate https://raw.githubusercontent.com/3dgiordano/SARS-CoV-2-Variants/main/data/variants.csv -O SARSCoV2_variants_${DATE}.csv || \
	{ echo Using fallback from ./data/; \
	 rm SARSCoV2_variants_${DATE}.csv; \
	 cp !{workflow.projectDir}/data/variants_SARSCoV2/SARSCoV2_variants_fallback_*.csv .; }
	'''
}
