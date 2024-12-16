include { kraken2 } from './process/kraken2.nf' 
include { krona } from './process/krona.nf' 
include { download_database_kraken2 } from './process/download_database_kraken2.nf'
include { freyja; freyja_plot } from './process/freyja.nf'
include { lcs_plot; lcs_sc2; lcs_ucsc_markers_table } from './process/lcs_sc2'
include { hercules } from './process/hercules.nf'

workflow read_classification_wf {
    take:   
        fastq
    main: 

    // Check for human contamination
        // database download
        if (params.krakendb) { kraken_db = file("${params.krakendb}") }
        else  { download_database_kraken2(); kraken_db = download_database_kraken2.out } 

        // classification
        kraken2(fastq, kraken_db)

        // visuals
        krona(kraken2.out)
        
    emit:   
        kraken = kraken2.out
}

workflow read_screening_freyja_wf {
    take:
        alignment
    main:
        freyja(alignment)
        freya_result_ch = freyja.out.aggregate.map{it -> it[1]}.collectFile(name: 'freyja_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.lineagedir}/")
        freyja_plot(freya_result_ch)
    
    emit:
        freyja_results = freyja.out.aggregate
        freyja_plots = freyja_plot.out
}

workflow read_screening_lsc_wf {
    take:
        fastq
    main:
        // Metagenomic analysis
        // calculate mixed/ pooled samples using LCS, https://github.com/rvalieris/LCS
        if (params.lcs_variant_groups == 'default')     { lcs_variant_groups_ch = Channel.empty() } 
        else                                            { lcs_variant_groups_ch = Channel.fromPath("${params.lcs_variant_groups}", checkIfExists: true)}

        lcs_ucsc_markers_table(lcs_variant_groups_ch.ifEmpty([]))
        lcs_sc2(fastq.combine(lcs_ucsc_markers_table.out))

        lcs_result_ch = lcs_sc2.out.map{it -> it[1]}.collectFile(name: 'lcs_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.lineagedir}/")
        lcs_plot(lcs_result_ch, params.lcs_cutoff)

    emit:
        lcs_results = lcs_sc2.out
        lcs_plots = lcs_plot.out
}

workflow read_screening_hercules_wf {
    take:
        fastq
    main:
        hercules(fastq)
    //emit:
}