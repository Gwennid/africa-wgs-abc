### 20210520
### Gwenna Breton
### Goal: make heatmaps and plots with the votes for the various analyses comparing continuous and pulse migration by topology. I will make plot for the models & groups analyses.

setwd("/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/modelchoice_pulseVScontinuous/QC/Scenario_by_scenario/")

###
### Topology 1a.
###

##############
######6 models
data <- read.table(file="ABC-RF-1a/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

library(gplots)
library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-1a/Confusion_matrix_1a_6models.pdf",width=10,height=10,pointsize = 10)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Type of migration",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,cexCol = 2,cexRow = 2, srtCol=0
)
dev.off()
#srtCol=90 for column names perpendicular to the matrix
#adjCol = c(0.5,1)

### Plot the votes

datA2 <- c(316,295,290,0,11,8)  
datB2 <- c(307,359,251,1,11,71)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-1a/Votes_1a_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.655",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.6560667",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-1a/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-1a/Confusion_matrix_1a_6models_2groups.pdf",width=5,height=5,pointsize = 8)
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

datA2 <- c(871,129)  
datB2 <- c(881,119)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-1a/Votes_1a_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9696333",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9635",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 1b.
###

##############
######6 models
data <- read.table(file="ABC-RF-1b/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

library(gplots)
library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-1b/Confusion_matrix_1b_6models.pdf",width=10,height=10,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Type of migration",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=0,cexCol = 2,cexRow = 2
)
dev.off()
#adjCol = c(0.5,1),

### Plot the votes

datA2 <- c(563,243,125,2,4,63)  
datB2 <- c(542,291,115,2,4,46)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-1b/Votes_1b_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.69345",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.6471667",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-1b/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-1b/Confusion_matrix_1b_6models_2groups.pdf",width=5,height=5,pointsize = 8)
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

datA2 <- c(898,102)  
datB2 <- c(912,88)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-1b/Votes_1b_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9749667",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9771333",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 1c.
###

##############
######6 models
data <- read.table(file="ABC-RF-1c/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

library(gplots)
library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-1c/Confusion_matrix_1c_6models.pdf",width=10,height=10,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Type of migration",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=0,cexCol = 2,cexRow = 2
)
dev.off()

### Plot the votes

datA2 <- c(510,250,175,3,11,51)
datB2 <- c(538,253,153,2,5,49)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-1c/Votes_1c_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.6332667",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.6778667",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-1c/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-1c/Confusion_matrix_1c_6models_2groups.pdf",width=5,height=5,pointsize = 8)
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

datA2 <- c(898,102)  
datB2 <- c(918,82)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-1c/Votes_1c_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.99095",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9935667",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 2a.
###

##############
######6 models
data <- read.table(file="ABC-RF-2a/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

#library(gplots)
#library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-2a/Confusion_matrix_2a_6models.pdf",width=10,height=10,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 2a, 6 models",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=0,cexCol = 2,cexRow = 2
)
dev.off()

### Plot the votes

datA2 <- c(521,249,51,1,8,170)
datB2 <- c(508,299,47,2,5,139)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-2a/Votes_2a_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.7046167",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.6866167",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-2a/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-2a/Confusion_matrix_2a_6models_2groups.pdf",width=5,height=5,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 2a, 2 groups",xlab="RF-ABC predicted model",ylab="True model",
          margins =c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=360,adjCol = c(0.5,1),cexCol = 2,cexRow = 2
)
dev.off()


### Plot the votes

datA2 <- c(894,106)  
datB2 <- c(915,85)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-2a/Votes_2a_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9841167",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9838667",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 2b.
###

##############
######6 models
data <- read.table(file="ABC-RF-2b/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

# library(gplots)
# library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-2b/Confusion_matrix_2b_6models.pdf",width=10,height=10,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 2b, 6 models",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=0,cexCol = 2,cexRow = 2
)
dev.off()

### Plot the votes

datA2 <- c(697,165,36,1,4,97)
datB2 <- c(703,183,40,1,2,71)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-2b/Votes_2b_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.7051333",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.7025167",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-2b/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-2b/Confusion_matrix_2b_6models_2groups.pdf",width=5,height=5,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 2b, 2 groups",xlab="RF-ABC predicted model",ylab="True model",
          margins =c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=360,adjCol = c(0.5,1),cexCol = 2,cexRow = 2
)
dev.off()


### Plot the votes

datA2 <- c(894,106)  
datB2 <- c(900,100)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-2b/Votes_2b_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9792833",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9748667",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 2c.
###

##############
######6 models
data <- read.table(file="ABC-RF-2c/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

library(gplots)
library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-2c/Confusion_matrix_2c_6models.pdf",width=10,height=10,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 2c, 6 models",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=0,cexCol = 2,cexRow = 2
)
dev.off()

### Plot the votes

datA2 <- c(473,252,66,3,11,195)
datB2 <- c(490,279,57,1,10,163)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-2c/Votes_2c_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.67135",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.6510167",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-2c/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-2c/Confusion_matrix_2c_6models_2groups.pdf",width=5,height=5,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 2c, 2 groups",xlab="RF-ABC predicted model",ylab="True model",
          margins =c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=360,adjCol = c(0.5,1),cexCol = 2,cexRow = 2
)
dev.off()


### Plot the votes

datA2 <- c(825,175)  
datB2 <- c(865,135)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-2c/Votes_2c_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9754833",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9875167",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()

###
### Topology 3a.
###

##############
######6 models
data <- read.table(file="ABC-RF-3a/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

library(gplots)
library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-3a/Confusion_matrix_3a_6models.pdf",width=10,height=10,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 3a, 6 models",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=0,cexCol = 2,cexRow = 2
)
dev.off()

### Plot the votes

datA2 <- c(472,265,75,3,10,175)
datB2 <- c(458,297,77,1,9,158)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-3a/Votes_3a_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.7022333",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.6515333",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-3a/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-3a/Confusion_matrix_3a_6models_2groups.pdf",width=5,height=5,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 3a, 2 groups",xlab="RF-ABC predicted model",ylab="True model",
          margins =c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=360,adjCol = c(0.5,1),cexCol = 2,cexRow = 2
)
dev.off()


### Plot the votes

datA2 <- c(862,138)  
datB2 <- c(874,126)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-3a/Votes_3a_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9724667",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.97245",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 3b.
###

##############
######6 models
data <- read.table(file="ABC-RF-3b/confusionmatrix_6Models.txt",header=TRUE)

m <- as.matrix(data)
rownames(m) <- c("P high","P inter","P low","C high","C inter","C low")
colnames(m) <- c("P high","P inter","P low","C high","C inter","C low")

# library(gplots)
# library(RColorBrewer)

my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,2000,length=100), # for red
               seq(2001,4000,length=100),  # for yellow
               seq(4001,6000,length=100)) # for green

par(mfrow=c(1,1))

#Comment: the plot is not perfect because the column labels overlap with the matrix. But it will do for now.
pdf(file="ABC-RF-3b/Confusion_matrix_3b_6models.pdf",width=10,height=10,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 3b, 6 models",xlab="RF-ABC predicted model",ylab="True model",
          margins=c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=0,cexCol = 2,cexRow = 2
)
dev.off()

### Plot the votes

datA2 <- c(650,202,38,1,7,102)
datB2 <- c(632,239,46,0,3,80)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-3b/Votes_3b_6models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2),mar=c(9,4,4,2)+0.1,oma=c(2,0,0,0))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.7195667",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-.5, "Observed: A2 (RHGn: Nzime)",outer=TRUE)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.7453333",ylab="Frequency",
        names.arg =c("pulse high","pulse inter","pulse low","continuous high","continuous inter","continuous low"),col=c(rep(col5[1],3),rep(col5[2],3)),cex.names=1,
        las=2)
mtext(side=1, line=8, at=0, adj=-.5, "Type of migration")
mtext(side=1, line=0, at=0, adj=-2, "Observed: B2 (RHGn: Ba.Konjo)",outer=TRUE)
dev.off()

##############
######2 groups

data <- read.table(file="ABC-RF-3b/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,16000,length=100), # for red
               seq(16001,32000,length=100),  # for yellow
               seq(32001,50000,length=100)) # for green

par(mfrow=c(1,1))

pdf(file="ABC-RF-3b/Confusion_matrix_3b_6models_2groups.pdf",width=5,height=5,pointsize = 8)
heatmap.2(m,
          cellnote = m,notecol="black",density.info="none",trace="none",Colv="NA",Rowv="NA",
          dendrogram = "none",key = FALSE,
          main="Confusion table topology 3b, 2 groups",xlab="RF-ABC predicted model",ylab="True model",
          margins =c(8,10),lwid=c(1,7),lhei=c(1,7),
          col=my_palette,breaks=col_breaks,
          notecex=2,srtCol=360,adjCol = c(0.5,1),cexCol = 2,cexRow = 2
)
dev.off()


### Plot the votes

datA2 <- c(903,97)  
datB2 <- c(915,85)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-3b/Votes_3b_6models_2groups.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9788333",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9827833",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()









setwd("/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/modelchoice_pulseVScontinuous/QC/Scenario_by_scenario/")



