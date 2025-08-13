process covarplot {
    label "covarplot"
    publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(vcf), path(depth1), path(depth2), path(primerbed)
    output:
        tuple val(name), path("${name}_amplicon_coverage.png"), path("${name}_amplicon_coverage_log.png")
    script:
        """
        covarplot.py -v ${vcf} -d1 ${depth1} -d2 ${depth2} -b ${primerbed}/${params.schemeLength == 400 ? 'artic' : 'varvamp'}-sars-cov-2/${params.schemeLength}/${params.primerV}/primer.bed -s .
        mv ${name}.CoVarPlot.png ${name}_amplicon_coverage.png
        covarplot.py -v ${vcf} -d1 ${depth1} -d2 ${depth2} -b ${primerbed}/${params.schemeLength == 400 ? 'artic' : 'varvamp'}-sars-cov-2/${params.schemeLength}/${params.primerV}/primer.bed -s . --log
        mv ${name}.CoVarPlot.png ${name}_amplicon_coverage_log.png
        """
    stub:
        """
        touch ${name}_amplicon_coverage.png ${name}_amplicon_coverage_log.png
        """
}

process covarplot_custom_bed {
    label "covarplot"
    publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy'
    input:
        tuple val(name), path(vcf), path(depth1), path(depth2), path(primerbed)
    output:
        tuple val(name), path("${name}_amplicon_coverage.png"), path("${name}_amplicon_coverage_log.png")
    script:
        """
        # clean up bed file: replace first colum with MN908947.3, remove empty lines and sort by 4th column (primer names) 
        cut -f2- ${primerbed} |\
            sed '/^[[:space:]]*\$/d' |\
            sed -e \$'s/^/MN908947.3\\t/' |\
            sort -k4 > nCoV-2019-plot.scheme.bed


        covarplot.py -v ${vcf} -d1 ${depth1} -d2 ${depth2} -b nCoV-2019-plot.scheme.bed -s .
        mv ${name}.CoVarPlot.png ${name}_amplicon_coverage.png
        covarplot.py -v ${vcf} -d1 ${depth1} -d2 ${depth2} -b nCoV-2019-plot.scheme.bed -s . --log
        mv ${name}.CoVarPlot.png ${name}_amplicon_coverage_log.png
        """
    stub:
        """
        touch ${name}_amplicon_coverage.png ${name}_amplicon_coverage_log.png
        """
}

/* USAGE
python3 CoVarPlot/covarplot.py \
-v work/0d/1e05c63f2b752e8b7ab75959623a4a/SAMPLE.pass.vcf.gz \
-d1 work/0d/1e05c63f2b752e8b7ab75959623a4a/SAMPLE.coverage_mask.txt.nCoV-2019_1.depths \
-d2 work/0d/1e05c63f2b752e8b7ab75959623a4a/SAMPLE.coverage_mask.txt.nCoV-2019_2.depths \
-b /home/hoelzerm/git/poreCov/data/external_primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed
*/
