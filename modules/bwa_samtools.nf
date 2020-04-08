process bwa_samtools {
    label "bwa"
    input:
        tuple val(name), path(fasta), path(reads)
    output:
        tuple val(name), path("coverage_info.txt"), optional: true
    script:
        """
        bwa index ${fasta}
        bwa mem ${fasta} ${reads} | samtools view -bS - | samtools sort -@ ${task.cpus} - > ${name}.1.bam
        samtools index -@ ${task.cpus} ${name}.1.bam

        samtools mpileup ${name}.1.bam | awk '{print \$1"\\t"\$2"\\t"\$4}' > coverage_info.txt

        if [ ! -s coverage_info.txt ] ; then
            rm coverage_info.txt
        fi
        """
}

