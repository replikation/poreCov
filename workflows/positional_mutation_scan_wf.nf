//include { positional_mutation_scan } from './process/positional_mutation_scan.nf'

workflow positional_mutation_scan_wf {
    take:   
        analysis_ch
        mutation_list
    main: 
        positional_mutation_scan(analysis_ch, mutation_list)
    emit:   
        positional_mutation_scan.out
}