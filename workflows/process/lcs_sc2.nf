process lcs_ucsc_markers_table {
    label 'lcs_sc2'

    input:
    path(variant_group_tsv)
    val(lsc_ucsc_work_version)

    output:
    path("LCS/outputs/variants_table/ucsc-markers-table-*.tsv")

    script:
    if ( params.lcs_ucsc_update || params.lcs_ucsc_version != 'predefined')
        """
        git clone https://github.com/rki-mf1/LCS.git --branch 2023.01.30

        if [[ "${params.lcs_variant_groups}" != default ]]; then
            rm -rf LCS/data/variant_groups.tsv
            cp ${variant_group_tsv} LCS/data/variant_groups.tsv
        fi

        cd LCS
        ## change settings
        sed -i "s/PB_VERSION=.*/PB_VERSION='${lsc_ucsc_work_version}'/" rules/config.py
        sed -i "s/NUM_SAMPLE=.*/NUM_SAMPLE=${params.lcs_ucsc_downsampling}/" rules/config.py
        mem=\$(echo ${task.memory} | cut -d' ' -f1)
        echo \$mem
        ## run pipeline
        snakemake --cores ${task.cpus} --resources mem_gb=\$mem --config dataset=somestring markers=ucsc -- ucsc_gather_tables
        ## output
        mv outputs/variants_table/ucsc-markers-table.tsv outputs/variants_table/ucsc-markers-table-${lsc_ucsc_work_version}.tsv 
        """
    else if ( params.lcs_ucsc_version == 'predefined' )
        """
        git clone https://github.com/rki-mf1/LCS.git --branch 2023.01.30
        mkdir -p LCS/outputs/variants_table
        zcat LCS/data/pre-generated-marker-tables/ucsc-markers-${params.lcs_ucsc_predefined}.tsv.gz > LCS/outputs/variants_table/ucsc-markers-table.tsv
        mv LCS/outputs/variants_table/ucsc-markers-table.tsv LCS/outputs/variants_table/ucsc-markers-table-predefined.tsv 
        """
    stub:
    """
    mkdir -p LCS/outputs/variants_table/
    touch LCS/outputs/variants_table/ucsc-markers-table-42.tsv
    """
}

process lcs_sc2 {
    label 'lcs_sc2'
    publishDir "${params.output}/${params.lineagedir}/${name}/lineage-proportion-by-reads", mode: 'copy'
    input:
    tuple val(name), path(reads), path(ucsc_markers_table)
  	output:
    tuple val(name), path("${name}.lcs.tsv")
  	script:
    """
    git clone https://github.com/rki-mf1/LCS.git --branch 2023.01.30

    mkdir -p LCS/outputs/variants_table
    mv ${ucsc_markers_table} LCS/outputs/variants_table/ucsc-markers-table.tsv  

    mkdir -p LCS/data/fastq
    cp ${reads} LCS/data/fastq/

    BN=\$(basename ${reads} .fastq.gz)
    echo \$BN >> LCS/data/tags_pool_mypool

    cd LCS
    snakemake --config markers=ucsc dataset=mypool --cores ${task.cpus} --resources mem_gb=8000 --set-threads pool_mutect=${task.cpus}
    cd ..

    cp LCS/outputs/decompose/mypool.out ${name}.lcs.tsv
    rm -rf LCS/data/fastq
    """
    stub:
    """
    touch ${name}.lcs.tsv
    """
  }

process lcs_plot {
  label "ggplot2"
  publishDir "${params.output}/${params.lineagedir}/", mode: 'copy'

  input:
  path(tsv)
  val(cutoff)
  
  output:
  path("*.png")
  
  script:
  """
  lcs_bar_plot.R '${tsv}' ${cutoff}
  """
  stub:
  """
  touch stub.png
  """
}