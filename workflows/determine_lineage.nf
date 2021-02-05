include { pangolin } from './process/pangolin'


workflow determine_lineage_wf {
    take: 
        fasta  
    main:
        pangolin(fasta)

        // collect lineage also to a summary     
        channel_tmp = pangolin.out.map {it -> it[1]}
                .splitCsv(header: true, sep: ',')
                .collectFile(seed: 'sequence_name,lineage,probability,pangoLEARN_version,status,note\n', 
                            storeDir: params.output + "/" + params.lineagedir + "/") {
                            row -> [ "metadata.tsv", row.taxon + ',' + row.lineage + ',' + row.probability + ',' + 
                            row.'pangoLEARN_version' + ',' + row.status + ',' + row.note + '\n']
                            }
    emit:
        pangolin.out
} 
