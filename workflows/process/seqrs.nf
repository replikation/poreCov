process seqrs {
        label 'seqrs'
        publishDir "${params.output}/${params.seqrepair}/", mode: 'copy'
    input:
        tuple val(name), path(fasta), path(primerbed)
  	output:
    	tuple val(name), file("*${name}.tsv")
  	script:
        """
        seqrs --genomes ${fasta} --primerbed ${primerbed}/v3.0.0/nCoV-2019.bed --results V3-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/v4.0.0/nCoV-2019.primer.bed --results V4-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/v4.1.0/nCoV-2019.primer.bed --results V4-1-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/v5.0.0/nCoV-2019.primer.bed --results V5-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/v5.1.0/nCoV-2019.primer.bed --results V5-1-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/v5.3.2/nCoV-2019.primer.bed --results V5.3.2_400-primer-to-repair-Ns-for_${name}.tsv -a 400
        seqrs --genomes ${fasta} --primerbed ${primerbed}/V1200/nCoV-2019.bed --results V1200-primer-to-repair-Ns-for_${name}.tsv -a 1200
        seqrs --genomes ${fasta} --primerbed ${primerbed}/v5.2.0/nCoV-2019.primer.bed --results V5.2.0_1200-primer-to-repair-Ns-for_${name}.tsv -a 1200
        """
    stub:
        """
        touch 1${name}.tsv
        """
}
