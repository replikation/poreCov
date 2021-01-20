include { toytree } from './process/toytree'

workflow toytree_wf {
    take: 
        trees  
    main:
        toytree(trees)
    emit:
        toytree.out
} 