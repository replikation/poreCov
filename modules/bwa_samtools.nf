process bwa_samtools {
    label "bwa"
    publishDir "${params.output}/fasta/${name}/", mode: 'copy', pattern: '*.sorted.bam*'
    input:
        tuple val(name), path(fasta), path(reads)
    output:
        tuple val(name), path("coverage_info.txt"), optional: true
        tuple val(name), path("*.sorted.bam"), path("*.sorted.bam.bai")
    script:
        """
        bwa index ${fasta}
        bwa mem ${fasta} ${reads} | samtools view -bS - | samtools sort -@ ${task.cpus} - > ${name}.sorted.bam
        samtools index -@ ${task.cpus} ${name}.sorted.bam

        samtools mpileup ${name}.sorted.bam | awk '{print \$1"\\t"\$2"\\t"\$4}' > coverage_info.txt

        if [ ! -s coverage_info.txt ] ; then
            rm coverage_info.txt
        fi
        """
}

