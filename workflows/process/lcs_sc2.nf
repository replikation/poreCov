process lcs_sc2 {
        label 'lcs_sc2'
        publishDir "${params.output}/${params.readqcdir}/${name}/mixed_sample_check", mode: 'copy'
    input:
        tuple val(name), path(reads)
  	output:
    	tuple val(name), path("${name}.lcs.out")
  	script:
    """
    git clone https://github.com/rvalieris/LCS.git

    mkdir -p LCS/outputs/variants_table
    zcat LCS/data/pre-generated-marker-tables/pango-designation-markers-v1.2.124.tsv.gz > LCS/outputs/variants_table/pango-markers-table.tsv
    zcat LCS/data/pre-generated-marker-tables/ucsc-markers-2022-01-31.tsv.gz > LCS/outputs/variants_table/ucsc-markers-table.tsv

    mkdir -p LCS/data/fastq
    cp ${reads} LCS/data/fastq/

    BN=\$(basename ${reads} .fastq.gz)
    echo \$BN >> LCS/data/tags_pool_mypool

    cd LCS
    snakemake --config markers=ucsc dataset=mypool --cores ${task.cpus} --resources mem_gb=8000
    cd ..

    cp LCS/outputs/decompose/mypool.out ${name}.lcs.out
    """
    stub:
    """
    touch ${name}.lcs.out
    """
  }