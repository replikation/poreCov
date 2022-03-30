#!/usr/bin/env Rscript
# via docker run --rm -it -v $PWD:/input r-base
#setwd("/input")

library(ggplot2)
library(dplyr)

args <- commandArgs(TRUE)
proportion_cutoff <- eval( parse(text=args[1]) )[1] 

input <- read.csv("lcs_results.tsv", sep='\t')

filter(input, proportion>=proportion_cutoff) %>% ggplot(aes(x=variant_group, y=proportion)) + geom_col(position='dodge') + 
geom_errorbar(width=0.4,aes(ymin=proportion-std_error, ymax=proportion+std_error)) + facet_wrap(~sample) +
theme(axis.text.x=element_text(angle=90,vjust=.5)) + coord_flip() + ggsave("lcs_bar_plot.png")

# for ggplot2 versions > 3.3.0:
# ggplot aes flip x and y
# geom_errorbar height
# rm coord_flip()