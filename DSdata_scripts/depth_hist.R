#!/usr/bin/Rscript 

args=commandArgs(trailingOnly=TRUE)

pdf(file=args[2], width=12, height=8)

a=read.table(args[1], head=F)
plot(a[,2], a[,3], cex=0.5)
hist(a[,3], breaks=500, cex=0.5)

dev.off()

