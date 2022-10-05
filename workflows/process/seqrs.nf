process seqrs {
        label 'seqrs'
        publishDir "${params.output}/${params.seqrepair}/", mode: 'copy'
    input:
        tuple val(name), path(fasta), path(primerbed)
  	output:
    	tuple val(name), file("*${name}.tsv")
  	script:
        """
        seqrs --genomes ${fasta} --primerbed ${primerbed}/V3/nCoV-2019.bed --results V3-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/V4/nCoV-2019.primer.bed --results V4-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/V4.1/nCoV-2019.primer.bed --results V4-1-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/V5/nCoV-2019.primer.bed --results V5-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/V1200/nCoV-2019.bed --results V1200-primer-to-repair-Ns-for_${name}.tsv -a 1200
        """
    stub:
        """
        touch 1${name}.tsv
        """
}
