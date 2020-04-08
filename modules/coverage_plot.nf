process coverage_plot {
    label "bwa"
    publishDir "${params.output}/coverage_info/", mode: 'copy'
    input:
        tuple val(name), path(coverage_info)
    output:
        tuple val(name), path("${name}_coverage.pdf"), path("${name}_coverage.svg")
    script:
        """
        coverage_plot.R
        mv overview.pdf ${name}_coverage.pdf 
        mv overview.svg ${name}_coverage.svg 
        """
}

