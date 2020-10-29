#!/usr/bin/env Rscript
# via docker run --rm -it -v $PWD:/input r-base
#setwd("/input")

library(ggplot2)
library(dplyr)
library(scales)

input  <- read.delim("coverage_info.txt", header = F, sep = "\t")


#png("phylum.png", height = 1000, width = 800, units = "px", pointsize = 12 )
samplename <- head(input, length(1))

plot <- ggplot(input) +
        geom_line(aes(x=V2, y=V3, color="coverage")) +
        theme_classic() +
        labs(title=as.vector(input$V1[[1]]), 
        x="Genome position",
        y="Coverage") +
        scale_y_continuous(trans = log10_trans(),
          breaks = trans_breaks("log10", function(x) 10^x),
          labels = trans_format("log10", math_format(10^.x)))


svg("overview.svg", height = 3, width = 12)
plot
dev.off()

pdf("overview.pdf", height = 3, width = 12)
plot
dev.off()


