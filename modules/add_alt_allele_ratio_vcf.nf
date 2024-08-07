process add_alt_allele_ratio_vcf {
    label 'artic'
    publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(bam), path(bai), path(vcf), path(failed_vcf)
        path(external_scheme) // primer scheme dir as input
    output:
       tuple val(name), path("${name}_all-vars-with-aar.vcf"), emit: vcf
       tuple val(name), path("mixed_sites_stats.csv"), emit: stats
    script:
    primer_version_tag = params.primerV.toString().contains(".bed") ? 'V_custom' : "${params.primerV}"
    primer_dir = params.primerV.toString().contains(".bed") ? "${external_scheme}" : "${external_scheme}/nCoV-2019"
    """
    # reheader failed VCF and change FILTER
    echo '##FILTER=<ID=ARTICFAIL,Description="ARTIC filter failed">' > add-to-hdr.txt
    # -c and --rename-annots add a header with a default/wrong description
    bcftools annotate -h add-to-hdr.txt -c "FILTER/ARTICFAIL:=FILTER/PASS" ${failed_vcf} | sed '/##FILTER=<ID=ARTICFAIL,Description="All filters passed">/d' > tmp_failed_updated-filter.vcf

    # concat failed and passed VCF
    bcftools concat ${vcf} tmp_failed_updated-filter.vcf | bcftools sort -o tmp_merged.vcf

    # call medaka tools annotate for each pool and add the alternate allele ratio
    for pool in `cut -f5 ${primer_dir}/${primer_version_tag}/nCoV-2019.scheme.bed | sort | uniq`; do
        bcftools view -i 'INFO/Pool="'\$pool'"' tmp_merged.vcf -o tmp_\$pool.vcf
        medaka tools annotate --dpsp --pad 25 --RG \$pool tmp_\$pool.vcf ${primer_dir}/${primer_version_tag}/nCoV-2019.reference.fasta ${bam} tmp_annotated_\$pool.vcf
        convert_VCF_info_fields.py tmp_annotated_\$pool.vcf tmp_aar_\$pool.vcf
    done

    # merge pool VCF and sort VCF
    bcftools concat tmp_aar_*.vcf | bcftools sort -o ${name}_all-vars-with-aar.vcf

    # count mixed sites
    # thresholds for human geotyping: 0.35 <= x <= 0.65 
    NUM_MIXED_SITES=\$(bcftools view -H -i 'INFO/DP>${params.min_depth} & INFO/AF>=0.3 & INFO/AF<=0.8' ${name}_all-vars-with-aar.vcf | wc -l)
    echo sample,num_mixed_sites > mixed_sites_stats.csv
    echo ${name},\$NUM_MIXED_SITES >> mixed_sites_stats.csv

    # remove intermediate files
    rm tmp_*
    """
    stub:
	"""
	touch ${name}_all-vars-with-aar.vcf mixed_sites_stats.csv
	"""
    
}
