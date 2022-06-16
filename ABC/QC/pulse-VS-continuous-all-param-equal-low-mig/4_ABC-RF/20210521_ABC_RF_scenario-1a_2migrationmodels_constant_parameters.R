### 20210521
### Gwenna Breton
### Goal: Perform ABC-RF to compare migration with pulses VS continous. Most parameters fixed (topology, divergence times, population sizes) and migration rates are fixed (10^-5 for pulse and 10^-6 for continuous).
#Based on 20210430_ABC_RF_pulseVScontinuous_low.R
#Based on code from Paul Verdu (ABC_RF_Gwenna_PrelimRES_FINAL_vSynthesis_18Feb2020.R) and various iterations of my code (see in /home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/Second_trial/ABC_modelchoiceRF_paramestimation_code).

setwd("/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/modelchoice_pulseVScontinuous/QC/Constant_parameters/")

library(abcrf)
library(abc)

########### All Stats


nbSim <- 1000
nbModels <- 2

ReSmodindex <- read.table("1a/1000repeats/ABC_scenario_1a_2models_1000repeats.MODINDEX")
ReSsumstat <- read.table("1a/1000repeats/ABC_scenario_1a_2models_1000repeats.SUMSTATALL")
nbsumstat <- length(as.character(colnames(ReSsumstat)))
nbsumstat
modindex2 <- as.factor(ReSmodindex$x)

########### Pour les vraies data
#two sets of five populations and two ways of computing the statistics. A=with Nzime, B=with BaKonjo; 1=strategy 1, 2=strategy 2.

folder_observed="/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/admixture/2_Second_trial/DATA/obs_sumstats/"

TrueDataA2 <- read.table(paste(folder_observed,"/strategy2/5Ju-5Ka-5Nz-5Baka-5Ns_100windows/200312.5Ju-5Ka-5Nz-5Baka-5Ns.100windows.vectorof382sumstats",sep=""), header=TRUE)
TrueDataB2 <- read.table(paste(folder_observed,"/strategy2/5Ju-5Ka-5BaK-5Baka-5Ns_100windows/200312.5Ju-5Ka-5BaK-5Baka-5Ns.100windows.vectorof382sumstats",sep=""), header=TRUE)
TrueDataA1 <- read.table(paste(folder_observed,"/strategy1/5Ju-5Ka-5Nz-5Baka-5Ns_100windows/200312.5Ju-5Ka-5Nz-5Baka-5Ns.100windows.vectorof382sumstats",sep=""), header=TRUE)
TrueDataB1 <- read.table(paste(folder_observed,"/strategy1/5Ju-5Ka-5BaK-5Baka-5Ns_100windows/200312.5Ju-5Ka-5BaK-5Baka-5Ns.100windows.vectorof382sumstats",sep=""), header=TRUE)


# # Check order sumstats
# 
# ##A2##
# SimTEMPa <- as.character(colnames(ReSsumstat))
# TrueTEMPa <- as.character(colnames(TrueDataA2))
# TrueTEMPb <- as.numeric(TrueDataA2[1,])
# 
# Check1 <- SimTEMPa==TrueTEMPa
# which(Check1==FALSE)
# which(is.na(TrueTEMPb)==TRUE)
# # C'est OK ici !
# 
# ##B2##
# TrueTEMPa <- as.character(colnames(TrueDataB2))
# TrueTEMPb <- as.numeric(TrueDataB2[1,])
# 
# Check1 <- SimTEMPa==TrueTEMPa
# which(Check1==FALSE)
# which(is.na(TrueTEMPb)==TRUE)
# # C'est OK ici !
# 
# 
# ########### Identify issues in the sumstats from the simulations
# 
# for (i in 1:nbsumstat){
#   if (length(which(is.na(as.numeric(ReSsumstat[,i]))==TRUE)>0)) {
#     temp <- which(is.na(as.numeric(ReSsumstat[,i]))==TRUE)
#     write.table(temp, file=paste("LIST_SimNumberMISSING_SumStatNumber_",i,"_",as.character(colnames(ReSsumstat)[i]),"_500Sims2Mod5pops_constant_20210521.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
#   }
# 
# }
# #Apparently no issues!
# 
# ###########
# ########### PCA with two models
# ###########
# 
# #I decided to remove the SFS distance statistics from this plot. The order of removal of statistics is changed compared to previously so I am naming the object otherwise.
# 
# ReSsumstat_noSFSdist <-  ReSsumstat[,-c(98:142)]
# nbsumstat1 <- length(colnames(ReSsumstat_noSFSdist))
# TrueDataA2_noSFSdist <- TrueDataA2[,-c(98:142)]
# TrueDataB2_noSFSdist <- TrueDataB2[,-c(98:142)]
# TrueDataA1_noSFSdist <- TrueDataA1[,-c(98:142)]
# TrueDataB1_noSFSdist <- TrueDataB1[,-c(98:142)]
# 
# ReSsumstat.pca <- prcomp(ReSsumstat_noSFSdist,scale=TRUE,center=TRUE)
# 
# #Project the observed data
# TrueDataA1.sc <- scale(TrueDataA1_noSFSdist, center = ReSsumstat.pca$center)
# TrueDataA1.pred <- TrueDataA1.sc %*% ReSsumstat.pca$rotation
# TrueDataB1.sc <- scale(TrueDataB1_noSFSdist, center = ReSsumstat.pca$center)
# TrueDataB1.pred <- TrueDataB1.sc %*% ReSsumstat.pca$rotation
# TrueDataA2.sc <- scale(TrueDataA2_noSFSdist, center = ReSsumstat.pca$center)
# TrueDataA2.pred <- TrueDataA2.sc %*% ReSsumstat.pca$rotation
# TrueDataB2.sc <- scale(TrueDataB2_noSFSdist, center = ReSsumstat.pca$center)
# TrueDataB2.pred <- TrueDataB2.sc %*% ReSsumstat.pca$rotation
# 
# #Plot
# library(RColorBrewer)
# lotsof_colours = c(brewer.pal(3,"Blues"),brewer.pal(3,"Reds"))
# couleurs <- c(rep(lotsof_colours[3],500),rep(lotsof_colours[6],500))
# 
# var_explained <- ReSsumstat.pca$sdev^2/sum(ReSsumstat.pca$sdev^2)
# var_explained_pc <- round(var_explained,4)*100
# 
# par(mfrow=c(1,1))
# #PC1-2
# pdf("ABC-RF-constant-1a_1000repeats/PC1_vs_PC2_337sumstats_20210524_1a_2models.pdf",height=10,width=10)
# plot(ReSsumstat.pca$x[,1],ReSsumstat.pca$x[,2],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
#      ylab=paste("PC2: ",var_explained_pc[2],"%",sep=""),main="PC2 versus PC1")
# points(TrueDataA1.pred[1],TrueDataA1.pred[2],col="black",pch=8)
# points(TrueDataA2.pred[1],TrueDataA2.pred[2],col="black",pch=2)
# points(TrueDataB1.pred[1],TrueDataB1.pred[2],col="black",pch=5)
# points(TrueDataB2.pred[1],TrueDataB2.pred[2],col="black",pch=1)
# legend(legend=c("pulse","continuous",
#                 "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
#        pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
# dev.off()
# #Nice! The clouds seem to overlap :)
# 
# #PC1-3
# pdf("ABC-RF-constant-1a_1000repeats/PC1_vs_PC3_337sumstats_20210524_1a_2models.pdf",height=10,width=10)
# plot(ReSsumstat.pca$x[,1],ReSsumstat.pca$x[,3],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
#      ylab=paste("PC3: ",var_explained_pc[3],"%",sep=""),main="PC3 versus PC1")
# points(TrueDataA1.pred[1],TrueDataA1.pred[3],col="black",pch=8)
# points(TrueDataA2.pred[1],TrueDataA2.pred[3],col="black",pch=2)
# points(TrueDataB1.pred[1],TrueDataB1.pred[3],col="black",pch=5)
# points(TrueDataB2.pred[1],TrueDataB2.pred[3],col="black",pch=1)
# legend(legend=c("pulse","continuous",
#                 "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
#        pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
# dev.off()
# 
# #PC1-4
# pdf("ABC-RF-constant-1a_1000repeats/PC1_vs_PC4_337sumstats_20210524_1a_2models.pdf",height=10,width=10)
# plot(ReSsumstat.pca$x[,1],ReSsumstat.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
#      ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC1")
# points(TrueDataA1.pred[1],TrueDataA1.pred[4],col="black",pch=8)
# points(TrueDataA2.pred[1],TrueDataA2.pred[4],col="black",pch=2)
# points(TrueDataB1.pred[1],TrueDataB1.pred[4],col="black",pch=5)
# points(TrueDataB2.pred[1],TrueDataB2.pred[4],col="black",pch=1)
# legend(legend=c("pulse","continuous",
#                  "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
#         pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
# dev.off()
# 
# #PC2-3
# pdf("ABC-RF-constant-1a_1000repeats/PC2_vs_PC3_337sumstats_20210524_1a_2models.pdf",height=10,width=10)
# plot(ReSsumstat.pca$x[,2],ReSsumstat.pca$x[,3],pch=20,asp=1,col=couleurs,xlab=paste("PC2: ",var_explained_pc[2],"%",sep=""),
#      ylab=paste("PC3: ",var_explained_pc[3],"%",sep=""),main="PC3 versus PC2")
# points(TrueDataA1.pred[2],TrueDataA1.pred[3],col="black",pch=8)
# points(TrueDataA2.pred[2],TrueDataA2.pred[3],col="black",pch=2)
# points(TrueDataB1.pred[2],TrueDataB1.pred[3],col="black",pch=5)
# points(TrueDataB2.pred[2],TrueDataB2.pred[3],col="black",pch=1)
# legend(legend=c("pulse","continuous",
#                 "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
#        pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
# dev.off()
# 
# #PC2-4
# pdf("ABC-RF-constant-1a_1000repeats/PC2_vs_PC4_337sumstats_20210524_1a_2models.pdf",height=10,width=10)
# plot(ReSsumstat.pca$x[,2],ReSsumstat.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC2: ",var_explained_pc[2],"%",sep=""),
#      ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC2")
# points(TrueDataA1.pred[2],TrueDataA1.pred[4],col="black",pch=8)
# points(TrueDataA2.pred[2],TrueDataA2.pred[4],col="black",pch=2)
# points(TrueDataB1.pred[2],TrueDataB1.pred[4],col="black",pch=5)
# points(TrueDataB2.pred[2],TrueDataB2.pred[4],col="black",pch=1)
# legend(legend=c("pulse","continuous",
#                  "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
#         pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
# dev.off()
# 
# #PC3-4
# pdf("ABC-RF-constant-1a_1000repeats/PC3_vs_PC4_337sumstats_20210524_1a_2models.pdf",height=10,width=10)
# plot(ReSsumstat.pca$x[,3],ReSsumstat.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC3: ",var_explained_pc[3],"%",sep=""),
#      ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC3")
# points(TrueDataA1.pred[3],TrueDataA1.pred[4],col="black",pch=8)
# points(TrueDataA2.pred[3],TrueDataA2.pred[4],col="black",pch=2)
# points(TrueDataB1.pred[3],TrueDataB1.pred[4],col="black",pch=5)
# points(TrueDataB2.pred[3],TrueDataB2.pred[4],col="black",pch=1)
# legend(legend=c("pulse","continuous",
#                 "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
#        pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
# dev.off()

########### First ABC Check indirect for low variances:

# dataX <- data.frame(modindex2,ReSsumstat)
# 
# modelTest <- abcrf(modindex2~., data = dataX, ntree=500, paral = TRUE, ncores = 8 , sampsize=nbModels*nbSim)

#Error in lda.default(x, grouping, ...) :
#  variables 208 209 210 211 212 223 224 225 226 227 228 229 230 231 232 258 259 263 264 268 269 270 273 274 275 278 279 287 appear to be constant within groups
#A bunch more of variables are constant (13)! (258 259 263 264 268 269 270 273 274 275 278 279 287)

ReSsumstat1 <- ReSsumstat[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232,258,259,263,264,268,269,270,273,274,275,278,279,287)]

nbsumstat1 <- length(colnames(ReSsumstat1)) #354

########### Real data: remove the same 28 statistics

TrueDataA2_1 <- TrueDataA2[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232,258,259,263,264,268,269,270,273,274,275,278,279,287)]
length(colnames(TrueDataA2_1))== nbsumstat1

TrueDataB2_1 <- TrueDataB2[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232,258,259,263,264,268,269,270,273,274,275,278,279,287)]
#length(colnames(TrueDataB2_1))== nbsumstat1

TrueDataA1_1 <- TrueDataA1[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232,258,259,263,264,268,269,270,273,274,275,278,279,287)]
TrueDataB1_1 <- TrueDataB1[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232,258,259,263,264,268,269,270,273,274,275,278,279,287)]

###########
#Subset (removed the statistics based on distance between pairs of SNPs in the same SFS bin).
TrueDataA2_2 <- TrueDataA2_1[,-c(98:142)]
TrueDataB2_2 <- TrueDataB2_1[,-c(98:142)]
TrueDataA1_2 <- TrueDataA1_1[,-c(98:142)]
TrueDataB1_2 <- TrueDataB1_1[,-c(98:142)]
ReSsumstat2 <- ReSsumstat1[,-c(98:142)]



# ########### Goodness of fit
# # Here I do it on the dataset where the low variance statistics are excluded, but maybe I shouldn't (?).
# 
# GFITdataA2_2 <- gfit(target=TrueDataA2_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
# png(filename = "ABC-RF-constant-1a_1000repeats/goodness_of_fit_ReSsumstat2_TrueDataA2_2_100replicates_tol0.01_20210525_1a_2models.png"); plot(GFITdataA2_2); dev.off()
# A2_p <- summary(GFITdataA2_2)$pvalue
# A2_d <- summary(GFITdataA2_2)$dist.obs
# 
# GFITdataB2_2 <- gfit(target=TrueDataB2_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
# png(filename = "ABC-RF-constant-1a_1000repeats/goodness_of_fit_ReSsumstat1_TrueDataB2_2_100replicates_tol0.01_20210525_1a_2models.png"); plot(GFITdataB2_2); dev.off()
# B2_p <- summary(GFITdataB2_2)$pvalue #0.65
# B2_d <- summary(GFITdataB2_2)$dist.obs
# 
# 
# GFITdataA1_2 <- gfit(target=TrueDataA1_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
# png(filename = "ABC-RF-constant-1a_1000repeats/goodness_of_fit_ReSsumstat1_TrueDataA1_2_100replicates_tol0.01_20210525_1a_2models.png"); plot(GFITdataA1_2); dev.off()
# A1_p <- summary(GFITdataA1_2)$pvalue
# A1_d <- summary(GFITdataA1_2)$dist.obs
# 
# 
# GFITdataB1_2 <- gfit(target=TrueDataB1_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
# png(filename = "ABC-RF-constant-1a_1000repeats/goodness_of_fit_ReSsumstat1_TrueDataB1_2_100replicates_tol0.01_20210525_1a_2models.png"); plot(GFITdataB1_2); dev.off()
# B1_p <- summary(GFITdataB1_2)$pvalue
# B1_d <- summary(GFITdataB1_2)$dist.obs
# 
# 
# # Make a summary plot about these results, with the observed distance and the p-value (Cf MetHis Suppl Figure S5).
# svg(file="ABC-RF-constant-1a_1000repeats/goodness_of_fit_ReSsumstat2_A1-A2-B1-B2_100replicates_tol0.01_20210525_1a_2models.svg",height = 12, width = 12)
# par(mfrow=c(2,2))
# plot(GFITdataA1_2, main=paste("Goodness-of-fit Nzime strategy 1 = ",round(A1_d,2), " (p-value=",A1_p,")",sep="")) #See how this works. It would be good to make it look similar to the plots of the summary statistics distributions.
# plot(GFITdataA2_2,main=paste("Goodness-of-fit Nzime strategy 2: ",round(A2_d,2), " (p-value=",A2_p,")",sep=""))
# plot(GFITdataB1_2,main=paste("Goodness-of-fit Ba.Konjo strategy 1: ",round(B1_d,2), " (p-value=",B1_p,")",sep=""))
# plot(GFITdataB2_2,main=paste("Goodness-of-fit Ba.Konjo strategy 2: ",round(B2_d,2), " (p-value=",B2_p,")",sep=""))
# dev.off()
# #TODO modify it manually (Inkscape) to match the colours in the sumstats distribution plots.





###########
########### Test ABC RF without real data for ReSsumstat2
###########

dataX2<- data.frame(modindex2,ReSsumstat2)
# modelTest2 <- abcrf(modindex2~., data = dataX2, ntree=1000, lda=FALSE, paral = TRUE, ncores = 8 , sampsize=nbModels*nbSim)
nbsumstat2 = length(colnames(ReSsumstat2)) #309

# png(file="ABC-RF-constant-1a_1000repeats/sumstats_importance_309sumstats_6Mod.png"); plot(modelTest2, dataX2, n.var=nbsumstat2, cex=0.5, pch=1); dev.off()
# png(file="ABC-RF-constant-1a_1000repeats/sumstats_importance_309sumstats_6Mod_TOP50.png"); plot(modelTest2, dataX2, n.var=50, cex=0.5, pch=1); dev.off()
# 
# modelTest2
# Call:
#   abcrf(formula = modindex2 ~ ., data = dataX2, lda = FALSE, ntree = 1000, sampsize = nbModels * nbSim, paral = TRUE, ncores = 8) 
# Number of simulations: 2000
# Out-of-bag prior error rate: 0.75%
# 
# Confusion matrix:
#   57  65 class.error
# 57 999   1       0.001
# 65  14 986       0.014
#The error rate is still super low! I had not expected that. Weird.

#Out-of-bag error of the forest for different numbers of trees:
# png(filename="ABC-RF-constant-1a_1000repeats/error_of_the_forest.png"); err.abcrf(modelTest2, training=dataX2, paral = TRUE, ncores = 8); dev.off()
#That does not look very good.

###########
########### Test ABC RF with real data for ReSsumstat2
###########

################## Sans LDA scores

# predDataA2_2 <- predict(modelTest2, obs=TrueDataA2_2, training=dataX2, ntree=1000, paral = TRUE, ncores = 8, sampsize=nbModels*nbSim)
# predDataA2_2
# # selected model votes model1 votes model2 post.proba
# # 1             65          489          511    0.99125
# 
# predDataB2_2 <- predict(modelTest2, obs=TrueDataB2_2, training=dataX2, ntree=1000, paral = TRUE, ncores = 8, sampsize=nbModels*nbSim)
# predDataB2_2
# # selected model votes model1 votes model2 post.proba
# # 1             65          495          505  0.9908833



###########
########### PCA with two models and only the summary statistics used for the model choice procedure.
###########

ReSsumstat2.pca <- prcomp(ReSsumstat2,scale=TRUE,center=TRUE)

#Project the observed data
TrueDataA1_2.sc <- scale(TrueDataA1_2, center = ReSsumstat2.pca$center)
TrueDataA1_2.pred <- TrueDataA1_2.sc %*% ReSsumstat2.pca$rotation
TrueDataB1_2.sc <- scale(TrueDataB1_2, center = ReSsumstat2.pca$center)
TrueDataB1_2.pred <- TrueDataB1_2.sc %*% ReSsumstat2.pca$rotation
TrueDataA2_2.sc <- scale(TrueDataA2_2, center = ReSsumstat2.pca$center)
TrueDataA2_2.pred <- TrueDataA2_2.sc %*% ReSsumstat2.pca$rotation
TrueDataB2_2.sc <- scale(TrueDataB2_2, center = ReSsumstat2.pca$center)
TrueDataB2_2.pred <- TrueDataB2_2.sc %*% ReSsumstat2.pca$rotation

#Plot
library(RColorBrewer)
lotsof_colours = c(brewer.pal(3,"Blues"),brewer.pal(3,"Reds"))
couleurs <- c(rep(lotsof_colours[3],500),rep(lotsof_colours[6],500))

var_explained <- ReSsumstat2.pca$sdev^2/sum(ReSsumstat2.pca$sdev^2)
var_explained_pc <- round(var_explained,4)*100

par(mfrow=c(1,1))
#PC1-2
pdf("ABC-RF-constant-1a_1000repeats/PC1_vs_PC2_309sumstats_20210524_1a_2models.pdf",height=10,width=10)
plot(ReSsumstat2.pca$x[,1],ReSsumstat2.pca$x[,2],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
     ylab=paste("PC2: ",var_explained_pc[2],"%",sep=""),main="PC2 versus PC1")
points(TrueDataA1_2.pred[1],TrueDataA1_2.pred[2],col="black",pch=8)
points(TrueDataA2_2.pred[1],TrueDataA2_2.pred[2],col="black",pch=2)
points(TrueDataB1_2.pred[1],TrueDataB1_2.pred[2],col="black",pch=5)
points(TrueDataB2_2.pred[1],TrueDataB2_2.pred[2],col="black",pch=1)
legend(legend=c("pulse","continuous",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
dev.off()
#Very similar to when I have all the sumstats.

#PC1-3
pdf("ABC-RF-constant-1a_1000repeats/PC1_vs_PC3_309sumstats_20210524_1a_2models.pdf",height=10,width=10)
plot(ReSsumstat2.pca$x[,1],ReSsumstat2.pca$x[,3],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
     ylab=paste("PC3: ",var_explained_pc[3],"%",sep=""),main="PC3 versus PC1")
points(TrueDataA1_2.pred[1],TrueDataA1_2.pred[3],col="black",pch=8)
points(TrueDataA2_2.pred[1],TrueDataA2_2.pred[3],col="black",pch=2)
points(TrueDataB1_2.pred[1],TrueDataB1_2.pred[3],col="black",pch=5)
points(TrueDataB2_2.pred[1],TrueDataB2_2.pred[3],col="black",pch=1)
legend(legend=c("pulse","continuous",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC1-4
pdf("ABC-RF-constant-1a_1000repeats/PC1_vs_PC4_309sumstats_20210524_1a_2models.pdf",height=10,width=10)
plot(ReSsumstat2.pca$x[,1],ReSsumstat2.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
     ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC1")
points(TrueDataA1_2.pred[1],TrueDataA1_2.pred[4],col="black",pch=8)
points(TrueDataA2_2.pred[1],TrueDataA2_2.pred[4],col="black",pch=2)
points(TrueDataB1_2.pred[1],TrueDataB1_2.pred[4],col="black",pch=5)
points(TrueDataB2_2.pred[1],TrueDataB2_2.pred[4],col="black",pch=1)
legend(legend=c("pulse","continuous",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC2-3
pdf("ABC-RF-constant-1a_1000repeats/PC2_vs_PC3_309sumstats_20210524_1a_2models.pdf",height=10,width=10)
plot(ReSsumstat2.pca$x[,2],ReSsumstat2.pca$x[,3],pch=20,asp=1,col=couleurs,xlab=paste("PC2: ",var_explained_pc[2],"%",sep=""),
     ylab=paste("PC3: ",var_explained_pc[3],"%",sep=""),main="PC3 versus PC2")
points(TrueDataA1_2.pred[2],TrueDataA1_2.pred[3],col="black",pch=8)
points(TrueDataA2_2.pred[2],TrueDataA2_2.pred[3],col="black",pch=2)
points(TrueDataB1_2.pred[2],TrueDataB1_2.pred[3],col="black",pch=5)
points(TrueDataB2_2.pred[2],TrueDataB2_2.pred[3],col="black",pch=1)
legend(legend=c("pulse","continuous",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC2-4
pdf("ABC-RF-constant-1a_1000repeats/PC2_vs_PC4_309sumstats_20210524_1a_2models.pdf",height=10,width=10)
plot(ReSsumstat2.pca$x[,2],ReSsumstat2.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC2: ",var_explained_pc[2],"%",sep=""),
     ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC2")
points(TrueDataA1_2.pred[2],TrueDataA1_2.pred[4],col="black",pch=8)
points(TrueDataA2_2.pred[2],TrueDataA2_2.pred[4],col="black",pch=2)
points(TrueDataB1_2.pred[2],TrueDataB1_2.pred[4],col="black",pch=5)
points(TrueDataB2_2.pred[2],TrueDataB2_2.pred[4],col="black",pch=1)
legend(legend=c("pulse","continuous",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC3-4
pdf("ABC-RF-constant-1a_1000repeats/PC3_vs_PC4_309sumstats_20210524_1a_2models.pdf",height=10,width=10)
plot(ReSsumstat2.pca$x[,3],ReSsumstat2.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC3: ",var_explained_pc[3],"%",sep=""),
     ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC3")
points(TrueDataA1_2.pred[3],TrueDataA1_2.pred[4],col="black",pch=8)
points(TrueDataA2_2.pred[3],TrueDataA2_2.pred[4],col="black",pch=2)
points(TrueDataB1_2.pred[3],TrueDataB1_2.pred[4],col="black",pch=5)
points(TrueDataB2_2.pred[3],TrueDataB2_2.pred[4],col="black",pch=1)
legend(legend=c("pulse","continuous",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,2),8,2,5,1),col=c(lotsof_colours[3],lotsof_colours[6],rep("black",4)),cex=0.8,ncol=2)
dev.off()


#######
####### Plot the top variables for the model choice.
#######
pdf("ABC-RF-constant-1a_1000repeats/bial.tot_pulseVScontinuous_1a.pdf",height=10,width=10)
boxplot(ReSsumstat2$bial.tot ~ modindex2, main="Total number of biallelic sites", names=c("Pulse","Continuous"), ylab="Number of sites",xlab="Type of migration")
dev.off()

pdf("ABC-RF-constant-1a_1000repeats/wat.theta.WP_pulseVScontinuous_1a.pdf",height=10,width=10)
boxplot(ReSsumstat2$wat.theta.WP ~ modindex2, main="Watterson's theta in the wRHG", names=c("Pulse","Continuous"), ylab="Watterson's theta",xlab="Type of migration")
dev.off()

pdf("ABC-RF-constant-1a_1000repeats/wat.theta.EP_pulseVScontinuous_1a.pdf",height=10,width=10)
boxplot(ReSsumstat2$wat.theta.EP ~ modindex2, main="Watterson's theta in the eRHG", names=c("Pulse","Continuous"), ylab="Watterson's theta",xlab="Type of migration")
dev.off()

pdf("ABC-RF-constant-1a_1000repeats/wat.theta.SK_pulseVScontinuous_1a.pdf",height=10,width=10)
boxplot(ReSsumstat2$wat.theta.SK ~ modindex2, main="Watterson's theta in the sKS", names=c("Pulse","Continuous"), ylab="Watterson's theta",xlab="Type of migration")
dev.off()

pdf("ABC-RF-constant-1a_1000repeats/wat.theta.NK_pulseVScontinuous_1a.pdf",height=10,width=10)
boxplot(ReSsumstat2$wat.theta.NK ~ modindex2, main="Watterson's theta in the nKS", names=c("Pulse","Continuous"), ylab="Watterson's theta",xlab="Type of migration")
dev.off()

pdf("ABC-RF-constant-1a_1000repeats/wat.theta.NP_pulseVScontinuous_1a.pdf",height=10,width=10)
boxplot(ReSsumstat2$wat.theta.NP ~ modindex2, main="Watterson's theta in the RHGn", names=c("Pulse","Continuous"), ylab="Watterson's theta",xlab="Type of migration")
dev.off()


