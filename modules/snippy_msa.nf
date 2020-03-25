process snippy_msa {
    label "snippy"
    input:
        tuple val(name), path(sample), file(snippy_dirs)
    output:
        tuple val(name), file("clean.full.aln")
    script:
    def maxmemory = "${task.memory.toString().replaceAll(/[\sGB]/,'')}"
        """
        snippy-core --ref ${sample} --prefix results ${snippy_dirs} 
        
        snippy-clean_full_aln results.full.aln > clean.full.aln
        """
}