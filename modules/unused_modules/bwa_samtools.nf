process bwa_samtools {
    label "bwa"
    publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: '*.sorted.bam*'
    input:
        tuple val(name), path(fasta), path(reads)
    output:
        tuple val(name), path("coverage_info.txt"), optional: true
        tuple val(name), path("*.sorted.bam"), path("*.sorted.bam.bai")
    script:
        """
        bwa index ${fasta}
        bwa mem -t ${task.cpus} ${fasta} ${reads} | samtools view -bS - | samtools sort -@ ${task.cpus} - > ${name}.sorted.bam
        samtools index -@ ${task.cpus} ${name}.sorted.bam

        samtools mpileup ${name}.sorted.bam | awk '{print \$1"\\t"\$2"\\t"\$4}' > coverage_info.txt

        if [ ! -s coverage_info.txt ] ; then
            rm coverage_info.txt
        fi

        # report stats

        """
}

/*

# .coverage.tsv
samtools depth ${name}.sorted.bam > ${name}.coverage.tsv

# pipeline.version
echo ${workflow.version} > pipeline.version

# .bamstats.pipe.txt 
samtools flagstat {input.bam} | sed -e 's/ + /|/' | sed -e 's/ /|/' 1> ${name}.bamstats.pipe.txt 

# .fragsize.tsv
echo 0 1 > ${name}.fragsize.tsv
samtools view -F 4 ${name}.sorted.bam | cut -f9 1>> ${name}.fragsize.tsv





# pangolin raw data
* needs to be forked from pangolin  > ${name}.lineage.txt

*/


/*

rule createReport:
    input:
        coverage = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.coverage.tsv"), sample = SAMPLES.keys()),
        frag_size  = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.fragsize.tsv"), sample = SAMPLES.keys()),
        mapping_statistics = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.txt"), sample = SAMPLES.keys()),
        mapping_statistics_forR = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.pipe.txt"), sample = SAMPLES.keys()),
        version = os.path.join(PROJFOLDER, "pipeline.version")
    output:
        report = os.path.join(PROJFOLDER, "qc_report.html"),
        csv = os.path.join(DATAFOLDER["reporting"], "coverage_samples.csv")
    params:
        p_folder = PROJFOLDER,
        l_folder = DATAFOLDER["reporting"],
        run_id = REPORT_RUNID,
        tax_id = KRAKEN_TAX_ID,
        template = srcdir("../ncov_minipipe.Rmd")
    log:
        os.path.join(DATAFOLDER["logs"], "reporting", "reporting.log")
    conda:
        "../envs/r.yaml"
    threads:
        1
    shell:
        # maybe need to replace shell by r call
        r"""
            # create report
            VERSION=${workflow.version}
            Rscript -e "rmarkdown::render('ncov_minipipe.Rmd',
                                            params=list(proj_folder='{params.p_folder}', list_folder='{params.l_folder}', run_name='{params.run_id}', tax_id='{params.tax_id}', version='$VERSION'),
                                            output_file=file.path('qc_report.html'))" 
        """
*/