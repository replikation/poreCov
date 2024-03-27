process freyja {
        label 'freyja'
        publishDir "${params.output}/${params.lineagedir}/${name}/lineage-proportion-by-reads", mode: 'copy', pattern: "*"

    input:
        tuple val(name), path(bam_file), path(reference)
    output:
        tuple val(name), path("*_freyja_aggregate.tsv"), emit: aggregate
        tuple val(name), path("*_freyja_variants.tsv"), path("*_freyja_depths.tsv"), path("*_freyja_demix.tsv"), path("*_freyja_aggregate.tsv"), path("*_freyja_plot.svg"), path("*_freyja_plot.png"), emit: fullresults

    script:
        if      ( params.freyja_update )  { freyja_update_cmd = "freyja update" }
        else                                { freyja_update_cmd = " " }
        """
        sed -i "s/^>.*\$/>MN908947.3/" ${reference} #rename reference-sequence as Alignment is done with this sequence-name even so using the same reference freyja will detect deviating names and thro an error
        mkdir -p freyja_result/

        ${freyja_update_cmd}
        
        freyja variants ${bam_file} \
            --variants ${name}_freyja_variants.tsv \
            --depths ${name}_freyja_depths.tsv \
            --ref ${reference}

        freyja demix ${name}_freyja_variants.tsv \
            ${name}_freyja_depths.tsv \
            --output freyja_result/${name}_freyja_demix.tsv

        freyja aggregate freyja_result/ \
            --output ${name}_freyja_aggregate.tsv

        freyja plot ${name}_freyja_aggregate.tsv \
            --output ${name}_freyja_plot.svg \
            --lineages

        freyja plot ${name}_freyja_aggregate.tsv \
            --output ${name}_freyja_plot.png \
            --lineages

        mv freyja_result/${name}_freyja_demix.tsv \$PWD
        """
}

process freyja_plot {
        label 'freyja'
        publishDir "${params.output}/${params.lineagedir}/", mode: 'copy', pattern: "*"

    input:
        path(aggregate_file)
    output:
        tuple path("*_freyja_plot.svg"), path("*_freyja_plot.png"), emit: results

    script:
        """
        freyja plot ${aggregate_file} \
            --output combined_freyja_plot.svg \
            --lineages

        freyja plot ${aggregate_file} \
            --output combined_freyja_plot.png \
            --lineages
        """
}