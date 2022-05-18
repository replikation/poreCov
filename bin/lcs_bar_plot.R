#!/usr/bin/env Rscript
# via docker run --rm -it -v $PWD:/input r-base
#setwd("/input")

library(ggplot2)
library(dplyr)

args <- commandArgs(TRUE)
tsv <- args[1]
# name <- parse(text=args[2])
proportion_cutoff <- eval( parse(text=args[2]) )

input <- read.csv(tsv, sep='\t')

filter(input, proportion>=proportion_cutoff) %>% ggplot(aes(x=variant_group, y=proportion)) + geom_col(position='dodge') + 
geom_errorbar(width=0.4,aes(ymin=proportion-std_error, ymax=proportion+std_error)) + 
facet_grid(sample ~ ., scales = "free", space = "free") +
theme(strip.text.y = element_text(angle = 0)) +
# facet_wrap(~sample, scales = "free_y", ncol=1) +
theme(axis.text.x=element_text(angle=90,vjust=.5)) + coord_flip() 
ggsave(paste0('lcs_barplot', '.png'), height=14)

# for ggplot2 versions > 3.3.0:
# ggplot aes flip x and y
# geom_errorbar height
# rm coord_flip()