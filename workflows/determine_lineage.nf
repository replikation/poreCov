include { pangolin } from './process/pangolin' 


workflow determine_lineage_wf {
    take: 
        fasta  
    main:
        pangolin(fasta)

        // collect lineage also to a summary     
        channel_tmp = pangolin.out[1]
                .splitCsv(header: true, sep: ',')
                .collectFile(seed: 'taxon,lineage,probability,pangoLEARN_version,status,note\n', 
                            storeDir: params.output + "/summary/") {
                            row -> [ "metadata.tsv", row.taxon + ',' + row.lineage + ',' + row.probability + ',' + 
                            row.'pangoLEARN_version' + ',' + row.status + ',' + row.note + '\n']
                            }
    emit:
        pangolin.out[0]
} 