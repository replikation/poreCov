process filter_fastq_by_length {
        label 'ubuntu'
        publishDir "${params.output}/${params.readsdir}/filtered_reads/", mode: 'copy', pattern: "${name}_filtered.fastq.gz"
    input:
        tuple val(name), path(reads) 
    output:
	    tuple val(name), path("${name}_filtered.fastq.gz"), optional: true
    script:
    read_min_length = params.minLength
    read_max_length = params.maxLength

    if ( params.primerV.matches('V1200') || params.primerV.matches('V5.2.0_1200') ) {
        if ( !params.minLength ) { read_min_length = 100 }
        if ( !params.maxLength ) { read_max_length = 1500 }
    }
    else {
        if ( !params.minLength ) { read_min_length = 100 }
        if ( !params.maxLength ) { read_max_length = 700 }
    }
    
        // we skip cleanup if samples are provided as the join channel removes unused barcodes
        """
        case "${reads}" in
            *.fastq.gz ) 
                zcat ${reads} | awk 'NR%4==1{a=\$0} NR%4==2{b=\$0} NR%4==3{c=\$0} NR%4==0&&length(b)>=${read_min_length}{print a"\\n"b"\\n"c"\\n"\$0;}' | gzip > "${name}_filtered.fastq.gz"
            ;;
            *.fastq)
                zcat ${reads} | awk 'NR%4==1{a=\$0} NR%4==2{b=\$0} NR%4==3{c=\$0} NR%4==0&&length(b)>=${read_min_length}{print a"\\n"b"\\n"c"\\n"\$0;}' | gzip > "${name}_filtered.fastq.gz"
            ;;
        esac
        
        if [[ ${params.samples} == false ]] && [[ ! "${params.fastq}" == *".csv" ]]; then
            find . -name "${name}_filtered.fastq.gz" -type 'f' -size -1500k -delete
        fi
        """
    stub:
        """
        touch ${name}_filtered.fastq.gz
        """    
}
