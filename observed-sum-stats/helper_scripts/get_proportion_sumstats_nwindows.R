#!/usr/bin/Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)
pop <- args[1]
windows <- args[2]

data <- read.table(file=paste(pop,".",windows,".1.2Mbpwindows.bial.nomiss.renamed.newpos.sumstatsH",sep=""),header=TRUE) #There was a warning "incomplete final line".
bial_multi <- read.table(file=paste(pop,".",windows,".1.2Mbpwindows.bi.multi.tot.sumstatsH",sep=""),header=TRUE)
prop_nomiss = data$bial.tot.nomiss/bial_multi$bial.tot
prop_bial = data[2:6]/data$bial.tot.nomiss
prop_hom = data[7:11]/data$bial.tot.nomiss
prop_private = data[32:36]/data$bial.tot.nomiss
write.table(file=paste(pop,".",windows,".1.2Mbpwindows.prop.sumstatsH",sep=""),c(prop_nomiss,prop_bial,prop_hom,prop_private),col.names=c('prop.bial.nomiss','prop.bial.NK','prop.bial.SK','prop.bial.WP','prop.bial.EP','prop.bial.NP','prop.hom.NK','prop.hom.SK','prop.hom.WP','prop.hom.EP','prop.hom.NP','prop.private.bial.NK','prop.private.bial.SK','prop.private.bial.WP','prop.private.bial.EP','prop.private.bial.NP'),quote=FALSE,row.names=FALSE)
