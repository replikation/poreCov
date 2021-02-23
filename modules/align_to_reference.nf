
process align_to_reference {
  	container = 'nanozoo/minimap2:2.17--7066fef'
    input:
        tuple(val(name), path(fastq_file), path(fasta_reference))
	output:
	    tuple(path("${name}.bam"), path("${name}.bam.bai"))
	script:
	"""
    minimap2 -ax map-ont ${fasta_reference} ${fastq_file} | samtools view -hbF4 | samtools sort > ${name}.bam
    samtools index ${name}.bam
	"""
}


