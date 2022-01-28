process artic_medaka {
        label 'artic'
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "*.consensus.fasta"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}_mapped_*.primertrimmed.sorted.bam*"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}.trimmed.rg.sorted.bam"
        publishDir "${params.output}/${params.genomedir}/all_consensus_sequences/", mode: 'copy', pattern: "*.consensus.fasta"

    input:
        tuple val(name), path(reads), path(external_scheme)
    output:
        tuple val(name), path("*.consensus.fasta"), emit: fasta
        tuple val(name), path("${name}_mapped_*.primertrimmed.sorted.bam"), path("${name}_mapped_*.primertrimmed.sorted.bam.bai"), emit: reference_bam
        tuple val(name), path("SNP_${name}.pass.vcf"), emit: vcf
        tuple val(name), path("${name}.pass.vcf.gz"), path("${name}.coverage_mask.txt.*1.depths"), path("${name}.coverage_mask.txt.*2.depths"), emit: covarplot
        tuple val(name), path("${name}.trimmed.rg.sorted.bam"), emit: fullbam
    script:   
        """
        artic minion    --medaka \
                        --medaka-model ${params.medaka_model} \
                        --min-depth ${params.min_depth} \
                        --normalise 500 \
                        --threads ${task.cpus} \
                        --scheme-directory ${external_scheme} \
                        --read-file ${reads} \
                        nCoV-2019/${params.primerV} ${name}

        # generate depth files
        artic_make_depth_mask --depth ${params.min_depth} \
            --store-rg-depths ${external_scheme}/nCoV-2019/${params.primerV}/nCoV-2019.reference.fasta \
            ${name}.primertrimmed.rg.sorted.bam \
            ${name}.coverage_mask.txt

        zcat ${name}.pass.vcf.gz > SNP_${name}.pass.vcf

        sed -i "1s/.*/>${name}/" *.consensus.fasta

        # get reference FASTA ID to rename BAM
        REF=\$(samtools view -H ${name}.primertrimmed.rg.sorted.bam | awk 'BEGIN{FS="\\t"};{if(\$1=="@SQ"){print \$2}}' | sed 's/SN://g')
        mv ${name}.primertrimmed.rg.sorted.bam ${name}_mapped_\${REF}.primertrimmed.sorted.bam
        samtools index ${name}_mapped_\${REF}.primertrimmed.sorted.bam
        """
        stub:
        """
        touch genome.consensus.fasta \
            ${name}_mapped_1.primertrimmed.sorted.bam \
            ${name}_mapped_1.primertrimmed.sorted.bam.bai \
            SNP_${name}.pass.vcf \
            ${name}.pass.vcf.gz \
            ${name}.coverage_mask.txt.nCoV-2019_1.depths \
            ${name}.coverage_mask.txt.nCoV-2019_2.depths \
            ${name}.trimmed.rg.sorted.bam
        """
}

process artic_nanopolish {
        label 'artic'
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "*.consensus.fasta"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}_mapped_*.primertrimmed.sorted.bam*"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}.trimmed.rg.sorted.bam"        
        publishDir "${params.output}/${params.genomedir}/all_consensus_sequences/", mode: 'copy', pattern: "*.consensus.fasta"

    input:
        tuple val(name), path(reads), path(external_scheme), path(fast5_dir), path(txt_files)
    output:
        tuple val(name), path("*.consensus.fasta"), emit: fasta
        tuple val(name), path("${name}_mapped_*.primertrimmed.sorted.bam"), path("${name}_mapped_*.primertrimmed.sorted.bam.bai"), emit: reference_bam
        tuple val(name), path("SNP_${name}.pass.vcf"), emit: vcf
        tuple val(name), path("${name}.pass.vcf.gz"), path("${name}.coverage_mask.txt.*1.depths"), path("${name}.coverage_mask.txt.*2.depths"), emit: covarplot
        tuple val(name), path("${name}.trimmed.rg.sorted.bam"), emit: fullbam
    script:   
        """
        artic minion --minimap2 --normalise 500 \
            --threads ${task.cpus} \
            --scheme-directory ${external_scheme} \
            --read-file ${reads} \
            --fast5-directory ${fast5_dir} \
            --sequencing-summary sequencing_summary*.txt \
            nCoV-2019/${params.primerV} ${name}

        # generate depth files
        artic_make_depth_mask --depth ${params.min_depth} \
            --store-rg-depths ${external_scheme}/nCoV-2019/${params.primerV}/nCoV-2019.reference.fasta \
            ${name}.primertrimmed.rg.sorted.bam \
            ${name}.coverage_mask.txt

        zcat ${name}.pass.vcf.gz > SNP_${name}.pass.vcf

        sed -i "1s/.*/>${name}/" *.consensus.fasta

        # get reference FASTA ID to rename BAM
        REF=\$(samtools view -H ${name}.primertrimmed.rg.sorted.bam | awk 'BEGIN{FS="\\t"};{if(\$1=="@SQ"){print \$2}}' | sed 's/SN://g')
        mv ${name}.primertrimmed.rg.sorted.bam ${name}_mapped_\${REF}.primertrimmed.sorted.bam
        samtools index ${name}_mapped_\${REF}.primertrimmed.sorted.bam
        """
        stub:
        """
        touch genome.consensus.fasta \
            ${name}_mapped_1.primertrimmed.sorted.bam \
            ${name}_mapped_1.primertrimmed.sorted.bam.bai \
            SNP_${name}.pass.vcf \
            ${name}.pass.vcf.gz \
            ${name}.coverage_mask.txt.nCoV-2019_1.depths \
            ${name}.coverage_mask.txt.nCoV-2019_2.depths \
            ${name}.trimmed.rg.sorted.bam
        """
}



// --min-depth
