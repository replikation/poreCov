process artic {
        label 'artic'
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy'
        //errorStrategy 'ignore'
    input:
        tuple val(name), path(reads), path(external_scheme)
    output:
        tuple val(name), path("*.consensus.fasta"), emit: fasta
        tuple val(name), path("SNP_${name}.pass.vcf"), emit: vcf
    script:
        """
        artic minion --medaka --medaka-model ${params.medaka_model} --normalise 200 --threads ${task.cpus} --scheme-directory ${external_scheme} \
            --read-file ${reads} nCoV-2019/${params.primerV} ${name}
        zcat ${name}.pass.vcf.gz > SNP_${name}.pass.vcf

        sed -i "1s/.*/>${name}/" *.consensus.fasta
        """
}