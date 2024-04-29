### 20210503
### Gwenna Breton
### Goal: Perform ABC-RF to compare migration with pulses VS continous. Priors: 0-1 and 0-0.05 respectively. 16 models.
#Based on 20210430_ABC_RF_pulseVScontinuous_low.R
#Based on code from Paul Verdu (ABC_RF_Gwenna_PrelimRES_FINAL_vSynthesis_18Feb2020.R) and various iterations of my code (see in /home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/Second_trial/ABC_modelchoiceRF_paramestimation_code).

setwd("/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/modelchoice_pulseVScontinuous/QC/Scenario_by_scenario/")

library(abcrf)
library(abc)

########### All Stats


nbSim <- 5000 
nbModels <- 6

ReSmodindex <- read.table("1a_analysisready/ABC_scenario_1a_6models_5000repeats.MODINDEX")
ReSsumstat <- read.table("1a_analysisready/ABC_scenario_1a_6models_5000repeats.SUMSTATALL")
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


# Check order sumstats

##A2##
SimTEMPa <- as.character(colnames(ReSsumstat))
TrueTEMPa <- as.character(colnames(TrueDataA2))
TrueTEMPb <- as.numeric(TrueDataA2[1,])

Check1 <- SimTEMPa==TrueTEMPa
which(Check1==FALSE)
which(is.na(TrueTEMPb)==TRUE)
# C'est OK ici !

##B2##
TrueTEMPa <- as.character(colnames(TrueDataB2))
TrueTEMPb <- as.numeric(TrueDataB2[1,])

Check1 <- SimTEMPa==TrueTEMPa
which(Check1==FALSE)
which(is.na(TrueTEMPb)==TRUE)
# C'est OK ici !


########### Identify issues in the sumstats from the simulations

# for (i in 1:nbsumstat){
#   if (length(which(is.na(as.numeric(ReSsumstat[,i]))==TRUE)>0)) {
#     temp <- which(is.na(as.numeric(ReSsumstat[,i]))==TRUE)
#     write.table(temp, file=paste("LIST_SimNumberMISSING_SumStatNumber_",i,"_",as.character(colnames(ReSsumstat)[i]),"_5000Sims16Mod5pops_low_20210517.txt", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
#   }
# 
# }
# #Apparently no issues!

########### First ABC Check indirect for low variances:

dataX <- data.frame(modindex2,ReSsumstat)

modelTest <- abcrf(modindex2~., data = dataX, ntree=500, paral = TRUE, ncores = 8 , sampsize=nbModels*nbSim)

#Error in lda.default(x, grouping, ...) :
#  variables 208 209 210 211 212 223 224 225 226 227 228 229 230 231 232 appear to be constant within groups
#Same list as in e.g. pulseVScontinuous high.

ReSsumstat1 <- ReSsumstat[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232)]

nbsumstat1 <- length(colnames(ReSsumstat1)) #367

########### Real data: remove the same 15 statistics

TrueDataA2_1 <- TrueDataA2[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232)]
length(colnames(TrueDataA2_1))== nbsumstat1

TrueDataB2_1 <- TrueDataB2[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232)]
#length(colnames(TrueDataB2_1))== nbsumstat1

TrueDataA1_1 <- TrueDataA1[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232)]
TrueDataB1_1 <- TrueDataB1[,-c(208,209,210,211,212,223,224,225,226,227,228,229,230,231,232)]

###########
#Subset (removed the statistics based on distance between pairs of SNPs in the same SFS bin).
TrueDataA2_2 <- TrueDataA2_1[,-c(98:142)]
TrueDataB2_2 <- TrueDataB2_1[,-c(98:142)]
TrueDataA1_2 <- TrueDataA1_1[,-c(98:142)]
TrueDataB1_2 <- TrueDataB1_1[,-c(98:142)]
ReSsumstat2 <- ReSsumstat1[,-c(98:142)]


###########
########### Test ABC RF without real data for ReSsumstat2
###########

dataX2<- data.frame(modindex2,ReSsumstat2)
modelTest2 <- abcrf(modindex2~., data = dataX2, ntree=1000, lda=FALSE, paral = TRUE, ncores = 8 , sampsize=nbModels*nbSim)
nbsumstat2 = length(colnames(ReSsumstat2)) #322

png(file="ABC-RF-1a/sumstats_importance_322sumstats_6Mod.png"); plot(modelTest2, dataX2, n.var=nbsumstat2, cex=0.5, pch=1); dev.off()
png(file="ABC-RF-1a/sumstats_importance_322sumstats_6Mod_TOP50.png"); plot(modelTest2, dataX2, n.var=50, cex=0.5, pch=1, pdf=TRUE); dev.off()

modelTest2
# Call:
#   abcrf(formula = modindex2 ~ ., data = dataX2, lda = FALSE, ntree = 1000, sampsize = nbModels * nbSim, paral = TRUE, ncores = 8) 
# Number of simulations: 30000
# Out-of-bag prior error rate: 14.2033%
# 
# Confusion matrix:
#   1   17   25   33   41   49 class.error
# 1  4144  677  113    0   19   47      0.1712
# 17  303 4337  358    0    0    2      0.1326
# 25   27  252 4721    0    0    0      0.0558
# 33    0    0    0 4062  936    2      0.1876
# 41    2    0    0 1393 3538   67      0.2924
# 49   30    8    0    0   25 4937      0.0126
#There is hardly any confusion between pulse and continuous. For the continuous, the lowest rate (model 49) stands out. Interestingly, it is confused with pulse high (which makes sense because it is the best from the continuous modalities, and pulse high is the best of models.)

# Out-of-bag error of the forest:
modelTest2$prior.err #14.20333% ; random: 83.3

#Out-of-bag error of the forest for different numbers of trees:
png(filename="ABC-RF-1a/error_of_the_forest.png"); err.abcrf(modelTest2, training=dataX2, paral = TRUE, ncores = 8); dev.off()

###########
########### Test ABC RF with real data for ReSsumstat2
###########

################## Sans LDA scores

predDataA2_2 <- predict(modelTest2, obs=TrueDataA2_2, training=dataX2, ntree=1000, paral = TRUE, ncores = 8, sampsize=nbModels*nbSim)
predDataA2_2
# selected model votes model1 votes model2 votes model3 votes model4 votes model5 votes model6 post.proba
# 1              1          316          295          290            0           11           88      0.655

predDataB2_2 <- predict(modelTest2, obs=TrueDataB2_2, training=dataX2, ntree=1000, paral = TRUE, ncores = 8, sampsize=nbModels*nbSim)
predDataB2_2
# selected model votes model1 votes model2 votes model3 votes model4 votes model5 votes model6 post.proba
# 1             17          307          359          251            1           11           71  0.6560667
#Model 17 is pulse intermediate.


###########
########### Group analysis
###########

#################################### 2 GROUPs = pulse versus continuous

modelTestgroup <- abcrf(modindex2~., data = dataX2,
                     group=list(c("1","17","25"),c("33","41","49")),
                     ntree=1000, lda=FALSE, paral = TRUE, ncores = 8 , sampsize=nbModels*nbSim)

png(file="ABC-RF-1a/sumstats_importance_322sumstats_2Groups.png"); plot(modelTestgroup, dataX2, n.var=nbsumstat2, cex=0.5, pch=1); dev.off()
png(file="ABC-RF-1a/sumstats_importance_322sumstats_2Groups_TOP50.png"); plot(modelTestgroup, dataX2, n.var=50, cex=0.5, pch=1, pdf=TRUE); dev.off()

modelTestgroup
# Call:
#   abcrf(formula = modindex2 ~ ., data = dataX2, group = list(c("1", "17", "25"), c("33", "41", "49")), lda = FALSE, ntree = 1000, sampsize = nbModels * nbSim, paral = TRUE, ncores = 8) 
# Number of simulations: 30000
# Out-of-bag prior error rate: 0.4633%
# 
# Confusion matrix:
#   g1    g2 class.error
# g1 14892   108 0.007200000
# g2    31 14969 0.002066667

# Predict
predDataA2_2b <- predict(modelTestgroup, obs=TrueDataA2_2, training=dataX2, ntree=1000, paral = TRUE, ncores = 8, sampsize=nbModels*nbSim)
predDataA2_2b
# selected group votes group1 votes group2 post.proba
# 1             g1          871          129  0.9696333
predDataB2_2b <- predict(modelTestgroup, obs=TrueDataB2_2, training=dataX2, ntree=1000, paral = TRUE, ncores = 8 , sampsize=nbModels*nbSim)
predDataB2_2b
# selected group votes group1 votes group2 post.proba
# 1             g1          881          119     0.9635

###########
########### PCA with six models
###########

#I decided to remove the SFS distance statistics from this plot. The order of removal of statistics is changed compared to previously so I am naming the object otherwise.

ReSsumstat_noSFSdist <-  ReSsumstat[,-c(98:142)]
nbsumstat1 <- length(colnames(ReSsumstat_noSFSdist))
TrueDataA2_noSFSdist <- TrueDataA2[,-c(98:142)]
TrueDataB2_noSFSdist <- TrueDataB2[,-c(98:142)]
TrueDataA1_noSFSdist <- TrueDataA1[,-c(98:142)]
TrueDataB1_noSFSdist <- TrueDataB1[,-c(98:142)]

ReSsumstat.pca <- prcomp(ReSsumstat_noSFSdist,scale=TRUE,center=TRUE)

#Project the observed data
TrueDataA1.sc <- scale(TrueDataA1_noSFSdist, center = ReSsumstat.pca$center)
TrueDataA1.pred <- TrueDataA1.sc %*% ReSsumstat.pca$rotation
TrueDataB1.sc <- scale(TrueDataB1_noSFSdist, center = ReSsumstat.pca$center)
TrueDataB1.pred <- TrueDataB1.sc %*% ReSsumstat.pca$rotation
TrueDataA2.sc <- scale(TrueDataA2_noSFSdist, center = ReSsumstat.pca$center)
TrueDataA2.pred <- TrueDataA2.sc %*% ReSsumstat.pca$rotation
TrueDataB2.sc <- scale(TrueDataB2_noSFSdist, center = ReSsumstat.pca$center)
TrueDataB2.pred <- TrueDataB2.sc %*% ReSsumstat.pca$rotation

#Plot
#Comment: It looks good! But I am not sure as to how to save it (large files...).
library(RColorBrewer)
lotsof_colours = c(brewer.pal(3,"Blues"),brewer.pal(3,"Reds"))
couleurs <- c(rep(lotsof_colours[1],5000),rep(lotsof_colours[2],5000),rep(lotsof_colours[3],5000),rep(lotsof_colours[4],5000),
              rep(lotsof_colours[5],5000),rep(lotsof_colours[6],5000))

var_explained <- ReSsumstat.pca$sdev^2/sum(ReSsumstat.pca$sdev^2)
var_explained_pc <- round(var_explained,4)*100

par(mfrow=c(1,1))
#PC1-2
pdf("ABC-RF-1a/PC1_vs_PC2_337sumstats_20210518_1a_6models.pdf",height=10,width=10)
plot(ReSsumstat.pca$x[,1],ReSsumstat.pca$x[,2],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
     ylab=paste("PC2: ",var_explained_pc[2],"%",sep=""),main="PC2 versus PC1")
points(TrueDataA1.pred[1],TrueDataA1.pred[2],col="black",pch=8)
points(TrueDataA2.pred[1],TrueDataA2.pred[2],col="black",pch=2)
points(TrueDataB1.pred[1],TrueDataB1.pred[2],col="black",pch=5)
points(TrueDataB2.pred[1],TrueDataB2.pred[2],col="black",pch=1)
legend(legend=c("pulse_high","pulse_intermediate","pulse_low","continuous_high","continuous_intermediate","continuous_low",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,6),8,2,5,1),col=c(lotsof_colours[1:6],rep("black",4)),cex=0.8,ncol=2)
dev.off()
#Nice! We actually see the cloud of continuous getting closer when the migration rates decrease.

#PC1-3
pdf("ABC-RF-1a/PC1_vs_PC3_337sumstats_20210518_1a_6models.pdf",height=10,width=10)
plot(ReSsumstat.pca$x[,1],ReSsumstat.pca$x[,3],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
     ylab=paste("PC3: ",var_explained_pc[3],"%",sep=""),main="PC3 versus PC1")
points(TrueDataA1.pred[1],TrueDataA1.pred[3],col="black",pch=8)
points(TrueDataA2.pred[1],TrueDataA2.pred[3],col="black",pch=2)
points(TrueDataB1.pred[1],TrueDataB1.pred[3],col="black",pch=5)
points(TrueDataB2.pred[1],TrueDataB2.pred[3],col="black",pch=1)
legend(legend=c("pulse_high","pulse_intermediate","pulse_low","continuous_high","continuous_intermediate","continuous_low",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,6),8,2,5,1),col=c(lotsof_colours[1:6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC1-4
pdf("ABC-RF-1a/PC1_vs_PC4_337sumstats_20210518_1a_6models.pdf",height=10,width=10)
plot(ReSsumstat.pca$x[,1],ReSsumstat.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC1: ",var_explained_pc[1],"%",sep=""),
     ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC1")
points(TrueDataA1.pred[1],TrueDataA1.pred[4],col="black",pch=8)
points(TrueDataA2.pred[1],TrueDataA2.pred[4],col="black",pch=2)
points(TrueDataB1.pred[1],TrueDataB1.pred[4],col="black",pch=5)
points(TrueDataB2.pred[1],TrueDataB2.pred[4],col="black",pch=1)
legend(legend=c("pulse_high","pulse_intermediate","pulse_low","continuous_high","continuous_intermediate","continuous_low",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,6),8,2,5,1),col=c(lotsof_colours[1:6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC2-3
pdf("ABC-RF-1a/PC2_vs_PC3_337sumstats_20210518_1a_6models.pdf",height=10,width=10)
plot(ReSsumstat.pca$x[,2],ReSsumstat.pca$x[,3],pch=20,asp=1,col=couleurs,xlab=paste("PC2: ",var_explained_pc[2],"%",sep=""),
     ylab=paste("PC3: ",var_explained_pc[3],"%",sep=""),main="PC3 versus PC2")
points(TrueDataA1.pred[2],TrueDataA1.pred[3],col="black",pch=8)
points(TrueDataA2.pred[2],TrueDataA2.pred[3],col="black",pch=2)
points(TrueDataB1.pred[2],TrueDataB1.pred[3],col="black",pch=5)
points(TrueDataB2.pred[2],TrueDataB2.pred[3],col="black",pch=1)
legend(legend=c("pulse_high","pulse_intermediate","pulse_low","continuous_high","continuous_intermediate","continuous_low",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,6),8,2,5,1),col=c(lotsof_colours[1:6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC2-4
pdf("ABC-RF-1a/PC2_vs_PC4_337sumstats_20210518_1a_6models.pdf",height=10,width=10)
plot(ReSsumstat.pca$x[,2],ReSsumstat.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC2: ",var_explained_pc[2],"%",sep=""),
     ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC2")
points(TrueDataA1.pred[2],TrueDataA1.pred[4],col="black",pch=8)
points(TrueDataA2.pred[2],TrueDataA2.pred[4],col="black",pch=2)
points(TrueDataB1.pred[2],TrueDataB1.pred[4],col="black",pch=5)
points(TrueDataB2.pred[2],TrueDataB2.pred[4],col="black",pch=1)
legend(legend=c("pulse_high","pulse_intermediate","pulse_low","continuous_high","continuous_intermediate","continuous_low",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,6),8,2,5,1),col=c(lotsof_colours[1:6],rep("black",4)),cex=0.8,ncol=2)
dev.off()

#PC3-4
pdf("ABC-RF-1a/PC3_vs_PC4_337sumstats_20210518_1a_6models.pdf",height=10,width=10)
plot(ReSsumstat.pca$x[,3],ReSsumstat.pca$x[,4],pch=20,asp=1,col=couleurs,xlab=paste("PC3: ",var_explained_pc[3],"%",sep=""),
     ylab=paste("PC4: ",var_explained_pc[4],"%",sep=""),main="PC4 versus PC3")
points(TrueDataA1.pred[3],TrueDataA1.pred[4],col="black",pch=8)
points(TrueDataA2.pred[3],TrueDataA2.pred[4],col="black",pch=2)
points(TrueDataB1.pred[3],TrueDataB1.pred[4],col="black",pch=5)
points(TrueDataB2.pred[3],TrueDataB2.pred[4],col="black",pch=1)
legend(legend=c("pulse_high","pulse_intermediate","pulse_low","continuous_high","continuous_intermediate","continuous_low",
                "Observed A1","Observed A2","Observed B1","Observed B2"),"topleft",
       pch=c(rep(20,6),8,2,5,1),col=c(lotsof_colours[1:6],rep("black",4)),cex=0.8,ncol=2)
dev.off()


########### Goodness of fit
# Here I do it on the dataset where the low variance statistics are excluded, but maybe I shouldn't (?).

GFITdataA2_2 <- gfit(target=TrueDataA2_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
png(filename = "ABC-RF-1a/goodness_of_fit_ReSsumstat2_TrueDataA2_2_100replicates_tol0.01_20210518_1a_6models.png"); plot(GFITdataA2_2); dev.off()
A2_p <- summary(GFITdataA2_2)$pvalue #0.41
A2_d <- summary(GFITdataA2_2)$dist.obs #71.72632

GFITdataB2_2 <- gfit(target=TrueDataB2_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
png(filename = "ABC-RF-1a/goodness_of_fit_ReSsumstat1_TrueDataB2_2_100replicates_tol0.01_20210518_1a_6models.png"); plot(GFITdataB2_2); dev.off()
B2_p <- summary(GFITdataB2_2)$pvalue #0.65
B2_d <- summary(GFITdataB2_2)$dist.obs


GFITdataA1_2 <- gfit(target=TrueDataA1_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
png(filename = "ABC-RF-1a/goodness_of_fit_ReSsumstat1_TrueDataA1_2_100replicates_tol0.01_20210518_1a_6models.png"); plot(GFITdataA1_2); dev.off()
A1_p <- summary(GFITdataA1_2)$pvalue
A1_d <- summary(GFITdataA1_2)$dist.obs


GFITdataB1_2 <- gfit(target=TrueDataB1_2, sumstat=ReSsumstat2, nb.replicate=100, tol=0.01, statistic=mean, subset=NULL, trace=FALSE)
png(filename = "ABC-RF-1a/goodness_of_fit_ReSsumstat1_TrueDataB1_2_100replicates_tol0.01_20210518_1a_6models.png"); plot(GFITdataB1_2); dev.off()
B1_p <- summary(GFITdataB1_2)$pvalue
B1_d <- summary(GFITdataB1_2)$dist.obs


# Make a summary plot about these results, with the observed distance and the p-value (Cf MetHis Suppl Figure S5).
svg(file="ABC-RF-1a/goodness_of_fit_ReSsumstat2_A1-A2-B1-B2_100replicates_tol0.01_20210518_1a_6models.svg",height = 12, width = 12)
par(mfrow=c(2,2))
plot(GFITdataA1_2, main=paste("Goodness-of-fit Nzime strategy 1 = ",round(A1_d,2), " (p-value=",A1_p,")",sep="")) #See how this works. It would be good to make it look similar to the plots of the summary statistics distributions.
plot(GFITdataA2_2,main=paste("Goodness-of-fit Nzime strategy 2: ",round(A2_d,2), " (p-value=",A2_p,")",sep=""))
plot(GFITdataB1_2,main=paste("Goodness-of-fit Ba.Konjo strategy 1: ",round(B1_d,2), " (p-value=",B1_p,")",sep=""))
plot(GFITdataB2_2,main=paste("Goodness-of-fit Ba.Konjo strategy 2: ",round(B2_d,2), " (p-value=",B2_p,")",sep=""))
dev.off()
#TODO modify it manually (Inkscape) to match the colours in the sumstats distribution plots.






