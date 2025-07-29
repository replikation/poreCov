include { pangolin } from './process/pangolin' 

workflow determine_lineage_wf {
    take: 
        fasta
        pangolindocker
    main:
        pangolin(fasta, pangolindocker)
        
        // collect lineage also to a summary     
        pangolin.out.map {it -> it[1]}
                .splitCsv(header: true, sep: ',')
                .collectFile(seed: 'taxon,lineage,conflict,pangolin_version,pangoLEARN_version,pango_version,status,note\n', 
                            storeDir: params.output + "/" + params.lineagedir + "/") {
                            row -> [ "metadata.csv", row.taxon + ',' + row.lineage + ',' + row.conflict + ',' + row.'pangolin_version' + ',' +
                            row.'pangoLEARN_version' + ',' + row.'pango_version' + ',' + row.status + ',' + row.note + '\n']
                            }
    emit:
        pangolin.out
} 
