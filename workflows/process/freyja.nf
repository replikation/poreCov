process freyja {
        label 'freyja'
        errorStrategy 'ignore'
        maxRetries 1
        publishDir "${params.output}/${params.lineagedir}/${name}/lineage-proportion-by-reads", mode: 'copy', pattern: "*"

    input:
        tuple val(name), path(bam_file), path(reference)
    output:
        tuple val(name), path("*_freyja_aggregate.tsv"), emit: aggregate
        tuple val(name), path("*_freyja_variants.tsv"), path("*_freyja_depths.tsv"), path("*_freyja_demix.tsv"), path("*_freyja_aggregate.tsv"), emit: fullresults
        tuple val(name), path("*_freyja_plot.svg"), path("*_freyja_plot.png"), emit: plots, optional: true

    script:
        if ( params.freyja_update )  
            {   freyja_update_cmd = "mkdir data; freyja update --outdir data" 
                // here we will download the latest Freyja data again instead of using the files obtained when installing the conda env
                freyja_ref_data = "data"
            }
        else
            {   freyja_update_cmd = " "
                // this is the original path where inside of the container the Freyja data is stored
                freyja_ref_data = "/opt/conda/envs/freyja-env/lib/python3.10/site-packages/freyja/data/"
            }
        """
        sed -i "s/^>.*\$/>MN908947.3/" ${reference} #rename reference-sequence as Alignment is done with this sequence-name even so using the same reference freyja will detect deviating names and throw an error
        mkdir -p freyja_result/

        ${freyja_update_cmd}

        freyja variants ${bam_file} \
            --variants ${name}_freyja_variants.tsv \
            --depths ${name}_freyja_depths.tsv \
            --ref ${reference}

        freyja demix ${name}_freyja_variants.tsv \
            ${name}_freyja_depths.tsv \
            --output freyja_result/${name}_freyja_demix.tsv \
            --barcodes ${freyja_ref_data}/usher_barcodes.csv \
            --lineageyml ${freyja_ref_data}/lineages.yml

        freyja aggregate freyja_result/ \
            --output ${name}_freyja_aggregate.tsv \

        freyja plot ${name}_freyja_aggregate.tsv \
            --output ${name}_freyja_plot.svg \
            --lineages \
            --mincov 10.0 #default minCoverage: 60 \
            --lineageyml ${freyja_ref_data}/lineages.yml

        freyja plot ${name}_freyja_aggregate.tsv \
            --output ${name}_freyja_plot.png \
            --lineages \
            --mincov 10.0 #default minCoverage: 60 \
            --lineageyml ${freyja_ref_data}/lineages.yml

        mv freyja_result/${name}_freyja_demix.tsv \$PWD #freyha_demix.tsv is put into a folder before to give it to freyja aggregate
        """
    stub:
        """
        touch stub_freyja_variants.tsv \
            stub_freyja_depths.tsv \
            stub_freyja_demix.tsv \
            stub_freyja_aggregate.tsv \
            stub_freyja_plot.png \
            stub_freyja_plot.svg
        """
}

process freyja_plot {
        label 'freyja'
        publishDir "${params.output}/${params.lineagedir}/", mode: 'copy', pattern: "*"

    input:
        path(aggregate_file)
    output:
        tuple path("*_freyja_plot.svg"), path("*_freyja_plot.png"), emit: results, optional: true

    script:
        """
        freyja plot ${aggregate_file} \
            --output combined_freyja_plot.svg \
            --lineages \
            --mincov 10.0 #default minCoverage: 60

        freyja plot ${aggregate_file} \
            --output combined_freyja_plot.png \
            --lineages \
            --mincov 10.0 #default minCoverage: 60
        """
    stub:
        """
        touch stub_combined_freyja_plot.png \
            stub_combined_freyja_plot.svg
        """
}