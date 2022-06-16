### 20210525
### Gwenna Breton
### Goal: make heatmaps and plots with the votes for the various analyses comparing continuous and pulse migration by topology. Most parameters are identical, only the migration modalities differ.
# Example for topology 1a.

setwd("/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/modelchoice_pulseVScontinuous/QC/Constant_parameters/")
library(gplots)
library(RColorBrewer)

###
### Topology 1a. Two models.
###

data <- read.table(file="ABC-RF-constant-1a_1000repeats/confusionmatrix_2Groups.txt",header=TRUE) #Replace input file.
#The input file is created manually (see example) from the output of the following command: "modelTest2 <- abcrf(modindex2~., data = dataX2, ntree=1000, lda=FALSE, paral = TRUE, ncores = 8 , sampsize=nbModels*nbSim)"

m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))

pdf(file="ABC-RF-constant-1a_1000repeats/Confusion_matrix_1a_2models.pdf",width=5,height=5,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Type of migration",xlab="RF-ABC predicted model",ylab="True model",
          margins =c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=360,adjCol = c(0.5,1),cexCol = 2,cexRow = 2
)
dev.off()


### Plot the votes
#The votes are the output of the command "predDataA2_2 <- predict(modelTest2, obs=TrueDataA2_2, training=dataX2, ntree=1000, paral = TRUE, ncores = 8, sampsize=nbModels*nbSim)"

datA2 <- c(489,511) #Replace by appropriate votes
datB2 <- c(495,505) #Replace by appropriate votes

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-constant-1a_1000repeats/Votes_1a_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.99125",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9908833",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()






