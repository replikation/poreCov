process filter_fastq_by_length {
        label 'ubuntu'
        publishDir "${params.output}/${params.readsdir}/filtered_reads/", mode: 'copy', pattern: "${name}_filtered.fastq.gz"
    input:
        tuple val(name), path(reads) 
    output:
	    tuple val(name), path("${name}_filtered.fastq.gz") optional true
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
                zcat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= ${read_min_length}' |\
                    awk -F"\\t" 'length(\$2)  <= ${read_max_length}' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
            ;;
            *.fastq)
                cat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= ${read_min_length}' |\
                awk -F"\\t" 'length(\$2)  <= ${read_max_length}' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
            ;;
        esac
        
        if [ ${params.samples} == false ]; then
            find . -name "${name}_filtered.fastq.gz" -type 'f' -size -1500k -delete
        fi
        """
    stub:
        """
        touch ${name}_filtered.fastq.gz
        """    
}
