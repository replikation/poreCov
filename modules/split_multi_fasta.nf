process split_multi_fasta {
    label 'ubuntu'
    input:
        tuple val(name), file(fasta) 
    output:
        tuple val(name), file("*.fa") 
    script:
    """
        cat ${fasta} | tr "\\ /|" "_" > all_clean_header.fasta

        while read line; do
            if [[ \${line:0:1} == '>' ]]
            then
                outfile=\${line#>}.fa
                echo "\${line}" > \${outfile}
            else
                echo "\${line}" >> \${outfile}
            fi
        done < all_clean_header.fasta
    """
}