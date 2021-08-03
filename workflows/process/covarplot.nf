process covarplot {
    label "covarplot"
    publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy'
    input:
        //tuple val(name), path(vcf), path(depth1), path(depth2), path(primer_bed)
        tuple val(name), path(vcf), path(primer_bed)
    output:
        tuple val(name), path("${name}_amplicon_coverage*.png")
    script:
        """
        covarplot.py -v ${vcf} -b ${primer_bed}/nCoV-2019/${params.primerV}/nCoV-2019.scheme.bed -s .
        mv ${name}.CoVarPlot.png ${name}_amplicon_coverage.png
        covarplot.py -v ${vcf} -b ${primer_bed}/nCoV-2019/${params.primerV}/nCoV-2019.scheme.bed -s . --log
        mv ${name}.CoVarPlot.png ${name}_amplicon_coverage_log.png
        """
    stub:
        """
        touch ${name}_amplicon_coverage1.png
        """
}

/* USAGE
python3 CoVarPlot/covarplot.py \
-v work/0d/1e05c63f2b752e8b7ab75959623a4a/SAMPLE.pass.vcf.gz \
-d1 work/0d/1e05c63f2b752e8b7ab75959623a4a/SAMPLE.coverage_mask.txt.nCoV-2019_1.depths \
-d2 work/0d/1e05c63f2b752e8b7ab75959623a4a/SAMPLE.coverage_mask.txt.nCoV-2019_2.depths \
-b /home/hoelzerm/git/poreCov/data/external_primer_schemes/nCoV-2019/V3/nCoV-2019.scheme.bed
*/
