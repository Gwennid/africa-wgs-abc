#20191213
#Gwenna
#Goal: compute summary statistics based on the pairwise allele sharing distance matrix. Script adapted to run on rackham and allows to specify name of input on the CL.

args <- commandArgs(trailingOnly = TRUE)
setwd(args[2])
ASD <- read.table(file=paste(args[1],".asd.dist",sep=""),header=TRUE,row.names = 1)

#I know that I will always have the same structure.

##
##Create vectors of pairwise ASD value for each population
##
#For NK
asdpopNK <- c()
for (i in c(1:5)) {for (j in c(1:5)) {if (i > j) {asdpopNK <- c(asdpopNK,ASD[i,j])}}}

#For SK
asdpopSK <- c()
for (i in c(6:10)) {for (j in c(6:10)) {if (i > j) {asdpopSK <- c(asdpopSK,ASD[i,j])}}}

#For WP
asdpopWP <- c()
for (i in c(11:15)) {for (j in c(11:15)) {if (i > j) {asdpopWP <- c(asdpopWP,ASD[i,j])}}}

#For EP
asdpopEP <- c()
for (i in c(16:20)) {for (j in c(16:20)) {if (i > j) {asdpopEP <- c(asdpopEP,ASD[i,j])}}}

#For NP
asdpopNP <- c()
for (i in c(21:25)) {for (j in c(21:25)) {if (i > j) {asdpopNP <- c(asdpopNP,ASD[i,j])}}}

##
## asd between each pair of populations
##
asdpopNKSK <- c()
for (i in c(6:10)) {for (j in c(1:5)) {asdpopNKSK <- c(asdpopNKSK,ASD[i,j])}}

asdpopNKWP <- c()
for (i in c(11:15)) {for (j in c(1:5)) {asdpopNKWP <- c(asdpopNKWP,ASD[i,j])}}

asdpopNKEP <- c()
for (i in c(16:20)) {for (j in c(1:5)) {asdpopNKEP <- c(asdpopNKEP,ASD[i,j])}}

asdpopNKNP <- c()
for (i in c(21:25)) {for (j in c(1:5)) {asdpopNKNP <- c(asdpopNKNP,ASD[i,j])}}

asdpopSKWP <- c()
for (i in c(11:15)) {for (j in c(6:10)) {asdpopSKWP <- c(asdpopSKWP,ASD[i,j])}}

asdpopSKEP <- c()
for (i in c(16:20)) {for (j in c(6:10)) {asdpopSKEP <- c(asdpopSKEP,ASD[i,j])}}

asdpopSKNP <- c()
for (i in c(21:25)) {for (j in c(6:10)) {asdpopSKNP <- c(asdpopSKNP,ASD[i,j])}}

asdpopWPEP <- c()
for (i in c(16:20)) {for (j in c(11:15)) {asdpopWPEP <- c(asdpopWPEP,ASD[i,j])}}

asdpopWPNP <- c()
for (i in c(21:25)) {for (j in c(11:15)) {asdpopWPNP <- c(asdpopWPNP,ASD[i,j])}}

asdpopEPNP <- c()
for (i in c(21:25)) {for (j in c(16:20)) {asdpopEPNP <- c(asdpopEPNP,ASD[i,j])}}

#Statistics based on projection of ASD

##If I project the entire matrix
projected = cmdscale(ASD,k=5)

#moyenne et variance intrapop pour chaque dimension
mean_NK = c(); mean_SK = c(); mean_WP = c(); mean_EP = c(); mean_NP = c()
var_NK = c(); var_SK = c(); var_WP = c(); var_EP = c(); var_NP = c()
for (i in c(1:5)) {
  mean_NK = c(mean_NK,mean(projected[1:5,i]))
  mean_SK = c(mean_SK,mean(projected[6:10,i]))
  mean_WP = c(mean_WP,mean(projected[11:15,i]))
  mean_EP = c(mean_EP,mean(projected[16:20,i]))
  mean_NP = c(mean_NP,mean(projected[21:25,i]))
  var_NK = c(var_NK,var(projected[1:5,i]))
  var_SK = c(var_SK,var(projected[6:10,i]))
  var_WP = c(var_WP,var(projected[11:15,i]))
  var_EP = c(var_EP,var(projected[16:20,i]))
  var_NP = c(var_NP,var(projected[21:25,i]))
}

#What kind of interpop values can I include here? Simply mean and variance for all individuals in each pair of population?
#for the mean I can simply take the average of the values in the two 'mean' vectors.
mean_NK_SK = c(); mean_NK_WP = c(); mean_NK_EP = c(); mean_NK_NP = c(); mean_SK_WP = c(); mean_SK_EP = c(); mean_SK_NP = c(); mean_WP_EP = c();
mean_WP_NP = c(); mean_EP_NP = c()
var_NK_SK = c(); var_NK_WP = c(); var_NK_EP = c(); var_NK_NP = c(); var_SK_WP = c(); var_SK_EP = c(); var_SK_NP = c(); var_WP_EP = c();
var_WP_NP = c(); var_EP_NP = c()

for (i in c(1:5)) {
  mean_NK_SK = c(mean_NK_SK,mean(projected[1:10,i]))
  mean_NK_WP = c(mean_NK_WP,mean(projected[c(1:5,11:15),i]))
  mean_NK_EP = c(mean_NK_EP,mean(projected[c(1:5,16:20),i]))
  mean_NK_NP = c(mean_NK_NP,mean(projected[c(1:5,21:25),i]))
  mean_SK_WP = c(mean_SK_WP,mean(projected[c(6:15),i]))
  mean_SK_EP = c(mean_SK_EP,mean(projected[c(6:10,16:20),i]))
  mean_SK_NP = c(mean_SK_NP,mean(projected[c(6:10,21:25),i]))
  mean_WP_EP = c(mean_WP_EP,mean(projected[c(11:20),i]))
  mean_WP_NP = c(mean_WP_NP,mean(projected[c(11:15,21:25),i]))
  mean_EP_NP = c(mean_EP_NP,mean(projected[c(16:25),i]))
  var_NK_SK = c(var_NK_SK,var(projected[1:10,i]))
  var_NK_WP = c(var_NK_WP,var(projected[c(1:5,11:15),i]))
  var_NK_EP = c(var_NK_EP,var(projected[c(1:5,16:20),i]))
  var_NK_NP = c(var_NK_NP,var(projected[c(1:5,21:25),i]))
  var_SK_WP = c(var_SK_WP,var(projected[c(6:15),i]))
  var_SK_EP = c(var_SK_EP,var(projected[c(6:10,16:20),i]))
  var_SK_NP = c(var_SK_NP,var(projected[c(6:10,21:25),i]))
  var_WP_EP = c(var_WP_EP,var(projected[c(11:20),i]))
  var_WP_NP = c(var_WP_NP,var(projected[c(11:15,21:25),i]))
  var_EP_NP = c(var_EP_NP,var(projected[c(16:25),i]))
}

#write out
write.table(file=paste(args[1],".asd.dist.sumstats",sep=""),x=t(c(mean(asdpopNK),mean(asdpopSK),mean(asdpopWP),mean(asdpopEP),mean(asdpopNP),var(asdpopNK),var(asdpopSK),var(asdpopWP),var(asdpopEP),var(asdpopNP),
  mean(asdpopNKSK),mean(asdpopNKWP),mean(asdpopNKEP),mean(asdpopNKNP),mean(asdpopSKWP),mean(asdpopSKEP),mean(asdpopSKNP),mean(asdpopWPEP),mean(asdpopWPNP),mean(asdpopEPNP),
  var(asdpopNKSK),var(asdpopNKWP),var(asdpopNKEP),var(asdpopNKNP),var(asdpopSKWP),var(asdpopSKEP),var(asdpopSKNP),var(asdpopWPEP),var(asdpopWPNP),var(asdpopEPNP),
  mean_NK,mean_SK,mean_WP,mean_EP,mean_NP,var_NK,var_SK,var_WP,var_EP,var_NP,
  mean_NK_SK,mean_NK_WP,mean_NK_EP,mean_NK_NP,mean_SK_WP,mean_SK_EP,mean_SK_NP,mean_WP_EP,mean_WP_NP,mean_EP_NP,
  var_NK_SK,var_NK_WP,var_NK_EP,var_NK_NP,var_SK_WP,var_SK_EP,var_SK_NP,var_WP_EP,var_WP_NP,var_EP_NP)),
  row.names=FALSE,col.names = FALSE,sep = "\t")

