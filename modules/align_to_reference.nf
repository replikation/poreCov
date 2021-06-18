process align_to_reference {
    label 'minimap2'
    input:
        tuple(val(name), path(fastq_file), path(fasta_reference))
    output:
        tuple(path("${name}.bam"), path("${name}.bam.bai"))
    script:
    """
    set -o pipefail
	
    minimap2 -t ${task.cpus} -ax map-ont ${fasta_reference} ${fastq_file} | samtools view -hbF4 | samtools sort -@ ${task.cpus} > ${name}.bam

    samtools index ${name}.bam
    """
    stub:
	"""
	touch ${name}.bam ${name}.bam.bai
	"""
    
}
