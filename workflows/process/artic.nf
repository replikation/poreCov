process artic {
        label 'artic'
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "*.consensus.fasta"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}_mapped_*.primertrimmed.sorted.bam*"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}.trimmed.rg.sorted.bam"
        publishDir "${params.output}/${params.genomedir}/all_consensus_sequences/", mode: 'copy', pattern: "*.consensus.fasta"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}.primersitereport.txt"
        publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "SNP_${name}.pass.vcf"
        publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}.coverage_mask.txt"
        publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}.fail.vcf"


    input:
        tuple val(name), path(reads), path(external_scheme)
        val(normalise_threshold)
    output:
        tuple val(name), path("*.consensus.fasta"), emit: fasta
        tuple val(name), path("${name}_mapped_*.primertrimmed.sorted.bam"), path("${name}_mapped_*.primertrimmed.sorted.bam.bai"), emit: reference_bam
        tuple val(name), path("SNP_${name}.pass.vcf"), emit: vcf
        tuple val(name), path("${name}.pass.vcf.gz"), path("${name}.coverage_mask.txt.*1.depths"), path("${name}.coverage_mask.txt.*2.depths"), emit: covarplot
        tuple val(name), path("${name}.trimmed.rg.sorted.bam"), emit: fullbam
        tuple val(name), path("${name}.primersitereport.txt"), emit: primersitereport
        tuple val(name), path("${name}.coverage_mask.txt"), emit: coverage_mask
        tuple val(name), path("${name}.fail.vcf"), emit: vcf_fail

    script:   
        def normalise_arg = normalise_threshold ? "--normalise ${normalise_threshold}" : '--normalise 0' // why is the --normalise flag not part of the bash script ^^
        
        // validate primer scheme version format
        def fetched_version = ("${params.primerV}" =~ ~ /[\d.]+/)[0]
        def formatted_version_number = null

        if (fetched_version.count('.') == 0){
            formatted_version_number = "v${fetched_version}.0.0"
        }
        else if (fetched_version.count('.') == 1){
            formatted_version_number = "v${fetched_version}.0"
        }
        else if (fetched_version.count('.') == 2){
            formatted_version_number = "v${fetched_version}"
        }
        else{
            formatted_version_number = null
            exit 1, println "\033[0;33mInvalid scheme version format 'v${fetched_version}' provided, please provide a version in the format 'vX.X.X', e.g. v1.0.0\u001B[0m"
        }

        """
        # artic minion --normalise 200 --scheme-directory ~/primer_schemes --scheme-name artic-inrb-mpox --scheme-length 2500 --scheme-version v1.0.0 --read-file run_name_barcode03.fastq samplename
        

        #  tolle wurst: Invalid scheme version format, please provide a version in the format 'vX.X.X', e.g. v1.0.0

        artic minion    --min-depth ${params.min_depth} \
                        ${normalise_arg} \
                        --threads ${task.cpus} \
                        --scheme-directory ${external_scheme} \
                        --read-file ${reads} \
                        --scheme-name nCov-2019 \
                        --scheme-version ${formatted_version_number} \
                        --scheme-length ${params.schemeLength} \
                        --model-dir ${params.clair3_model_dir} \
                        --model ${params.clair3_model_name} \
                        ${name}

        echo 'artic minion ran successfully'

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
            ${name}.trimmed.rg.sorted.bam \
            ${name}.primersitereport.txt \
            ${name}.coverage_mask.txt \
            ${name}.fail.vcf 
        """
}

process artic_custom_bed {
        label 'artic'
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "*.consensus.fasta"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}_mapped_*.primertrimmed.sorted.bam*"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}.trimmed.rg.sorted.bam"
        publishDir "${params.output}/${params.genomedir}/all_consensus_sequences/", mode: 'copy', pattern: "*.consensus.fasta"
        publishDir "${params.output}/${params.genomedir}/${name}/", mode: 'copy', pattern: "${name}.primersitereport.txt"
        publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "SNP_${name}.pass.vcf"
        publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}.coverage_mask.txt"
        publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}.fail.vcf"

    input:
        tuple val(name), path(reads), path(external_scheme), path(primerBed)
        val(normalise_threshold)
    output:
        tuple val(name), path("*.consensus.fasta"), emit: fasta
        tuple val(name), path("${name}_mapped_*.primertrimmed.sorted.bam"), path("${name}_mapped_*.primertrimmed.sorted.bam.bai"), emit: reference_bam
        tuple val(name), path("SNP_${name}.pass.vcf"), emit: vcf
        tuple val(name), path("${name}.pass.vcf.gz"), path("${name}.coverage_mask.txt.*1.depths"), path("${name}.coverage_mask.txt.*2.depths"), emit: covarplot
        tuple val(name), path("${name}.trimmed.rg.sorted.bam"), emit: fullbam
        tuple val(name), path("${name}.primersitereport.txt"), emit: primersitereport
        tuple val(name), path("${name}.coverage_mask.txt"), emit: coverage_mask
        tuple val(name), path("${name}.fail.vcf"), emit: vcf_fail
        path ("primer_scheme/nCoV-2019/"), emit: primer_dir
    script:   
        def normalise_arg = normalise_threshold ? "--normalise ${normalise_threshold}" : '--normalise 0'
        """
        # create a new primer dir as input for artic
        mkdir -p primer_scheme/nCoV-2019
        cp -r ${external_scheme}/nCoV-2019/V_custom primer_scheme/nCoV-2019/

        # clean up bed file: replace first colum with MN908947.3, remove empty lines and sort by 4th column (primer names) 
        cut -f2- ${primerBed} |\
            sed '/^[[:space:]]*\$/d' |\
            sed -e \$'s/^/MN908947.3\\t/' |\
            sort -k4 > primer_scheme/nCoV-2019/V_custom/nCoV-2019.scheme.bed

        # start artic
        artic minion    --medaka \
                        --medaka-model ${params.medaka_model} \
                        --min-depth ${params.min_depth} \
                        ${normalise_arg} \
                        --threads ${task.cpus} \
                        --scheme-directory primer_scheme \
                        --read-file ${reads} \
                        nCoV-2019/V_custom ${name}

        # generate depth files
        artic_make_depth_mask --depth ${params.min_depth} \
            --store-rg-depths primer_scheme/nCoV-2019/V_custom/nCoV-2019.reference.fasta \
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
            ${name}.trimmed.rg.sorted.bam \
            ${name}.primersitereport.txt \
            ${name}.coverage_mask.txt \
            ${name}.fail.vcf
        """
}


// --min-depth
