process filter_fastq_by_length {
        label 'ubuntu'
    input:
        tuple val(name), path(reads) 
    output:
	    tuple val(name), path("${name}_filtered.fastq.gz") optional true
    script:
    if ( params.primerV.matches('V1200') )
    """
    case "${reads}" in
        *.fastq.gz ) 
            zcat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= 250' | sed 's/\\t/\\n/g' |\
                awk -F"\\t" 'length(\$2)  <= 1500' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
        ;;
        *.fastq)
            cat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= 1000' | sed 's/\\t/\\n/g' |\
            awk -F"\\t" 'length(\$2)  <= 1300' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
        ;;
    esac

    find . -name "${name}_filtered.fastq.gz" -type 'f' -size -1500k -delete
    """
    else
    """
    case "${reads}" in
        *.fastq.gz ) 
            zcat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= ${params.minLength}' | sed 's/\\t/\\n/g' |\
                awk -F"\\t" 'length(\$2)  <= ${params.maxLength}' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
        ;;
        *.fastq)
            cat ${reads} | paste - - - - | awk -F"\\t" 'length(\$2)  >= ${params.minLength}' | sed 's/\\t/\\n/g' |\
            awk -F"\\t" 'length(\$2)  <= ${params.maxLength}' | sed 's/\\t/\\n/g' | gzip > "${name}_filtered.fastq.gz"
        ;;
    esac

    find . -name "${name}_filtered.fastq.gz" -type 'f' -size -1500k -delete
    """
}
