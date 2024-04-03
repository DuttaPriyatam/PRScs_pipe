#!/usr/bin/env Rscript

args= commandArgs(trailingOnly=T)
path <- args[1]
input <- paste(path,"raw.QC.het", sep='')
dat <- read.table(input, header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], paste(path,"raw.valid.sample",sep=''), quote=F, row.names=F) # print FID and IID for valid samples
