process align_to_reference {
    publishDir "${params.output}/bams/", mode: 'copy'
    label 'minimap2'
    input:
        tuple(val(name), path(fastq_file), path(fasta_reference))
    output:
        tuple(path("${name}.bam"), path("${name}.bam.bai"))
    script:
    """
	set -o pipefail
	
    minimap2 -ax map-ont ${fasta_reference} ${fastq_file} | samtools view -hbF4 | samtools sort > ${name}.bam

    samtools index ${name}.bam
    """
}