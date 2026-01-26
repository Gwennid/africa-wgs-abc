#!/usr/bin/Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)
pop <- args[1]
bi <- read.table(file=paste(file=pop,"_bi",sep=""))
multi <- read.table(file=paste(file=pop,"_multi",sep=""))
BI <- sum(bi)
MULTI <- sum(multi)
CONTIG <- nrow(bi)/4

write.table(file=paste(file=pop,".",CONTIG,".1.2Mbpwindows.bi.multi.tot.sumstats",sep=""),t(c(BI,MULTI)),row.names=FALSE,col.names=FALSE,sep="\t")

