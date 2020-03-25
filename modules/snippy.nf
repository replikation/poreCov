process snippy {
    label "snippy"
    input:
        tuple val(name), path(sample), val(reference_name), path(references)
    output:
        tuple val(name), path(sample), path("results_calling/${reference_name}")
    script:
    def maxmemory = "${task.memory.toString().replaceAll(/[\sGB]/,'')}"
        """
        snippy --cpus ${task.cpus} --ram ${maxmemory} --outdir results_calling/${reference_name} --ref ${sample} --ctgs ${references}
        """
}

