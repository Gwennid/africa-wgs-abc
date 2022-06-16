### 20210525
### Gwenna Breton
### Goal: make heatmaps and plots with the votes for the various analyses comparing continuous and pulse migration by topology. Most parameters are identical, only the migration modalities differ.

setwd("/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/modelchoice_pulseVScontinuous/QC/Constant_parameters/")
library(gplots)
library(RColorBrewer)

###
### Topology 1a. Two models.
###

data <- read.table(file="ABC-RF-constant-1a_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
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

datA2 <- c(489,511)  
datB2 <- c(495,505)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-constant-1a_1000repeats/Votes_1a_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.99125",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9908833",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 1b. Two models.
###

data <- read.table(file="ABC-RF-constant-1b_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))

pdf(file="ABC-RF-constant-1b_1000repeats/Confusion_matrix_1b_2models.pdf",width=5,height=5,pointsize = 8)
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

datA2 <- c(483,517)  
datB2 <- c(483,517)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-constant-1b_1000repeats/Votes_1b_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9921333",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9948",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 1c. Two models.
###

data <- read.table(file="ABC-RF-constant-1c_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))

pdf(file="ABC-RF-constant-1c_1000repeats/Confusion_matrix_1c_2models.pdf",width=5,height=5,pointsize = 8)
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

datA2 <- c(488,512)  
datB2 <- c(487,513)

col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.

pdf(file="ABC-RF-constant-1c_1000repeats/Votes_1c_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9828167",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.99275",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 2a. Two models.
###

data <- read.table(file="ABC-RF-constant-2a_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))

pdf(file="ABC-RF-constant-2a_1000repeats/Confusion_matrix_2a_2models.pdf",width=5,height=5,pointsize = 8)
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
datA2 <- c(527,473)
datB2 <- c(528,472)
col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.
pdf(file="ABC-RF-constant-2a_1000repeats/Votes_2a_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.98",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9885",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 2b. Two models.
###
data <- read.table(file="ABC-RF-constant-2b_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))
pdf(file="ABC-RF-constant-2b_1000repeats/Confusion_matrix_2b_2models.pdf",width=5,height=5,pointsize = 8)
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
datA2 <- c(483,517)
datB2 <- c(503,497)
col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.
pdf(file="ABC-RF-constant-2b_1000repeats/Votes_2b_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9811167",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9938",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 2c. Two models.
###
data <- read.table(file="ABC-RF-constant-2c_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))
pdf(file="ABC-RF-constant-2c_1000repeats/Confusion_matrix_2c_2models.pdf",width=5,height=5,pointsize = 8)
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
datA2 <- c(546,454)
datB2 <- c(554,446)
col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.
pdf(file="ABC-RF-constant-2c_1000repeats/Votes_2c_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.9928",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9963",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 3a. Two models.
###
data <- read.table(file="ABC-RF-constant-3a_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))
pdf(file="ABC-RF-constant-3a_1000repeats/Confusion_matrix_3a_2models.pdf",width=5,height=5,pointsize = 8)
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
datA2 <- c(520,480)
datB2 <- c(516,484)
col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.
pdf(file="ABC-RF-constant-3a_1000repeats/Votes_3a_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.98",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.98495",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()


###
### Topology 3b. Two models.
###
data <- read.table(file="ABC-RF-constant-3b_1000repeats/confusionmatrix_2Groups.txt",header=TRUE)
m <- as.matrix(data)
rownames(m) <- c("pulse","continuous")
colnames(m) <- c("pulse","continuous")
my_palette <- colorRampPalette(c("beige", "blue", "red"))(n = 299)
col_breaks = c(seq(0,400,length=100),
               seq(401,800,length=100),
               seq(801,1200,length=100))

par(mfrow=c(1,1))
pdf(file="ABC-RF-constant-3b_1000repeats/Confusion_matrix_3b_2models.pdf",width=5,height=5,pointsize = 8)
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
datA2 <- c(561,439)
datB2 <- c(533,467)
col5 <- brewer.pal(5, "Set2") #I will use the two first colours for the internal split plots, and the three last for the first branch to diverge.
pdf(file="ABC-RF-constant-3b_1000repeats/Votes_3b_2models.pdf",width=5,height=5,pointsize=8)
par(mfrow=c(1,2))
barplot(as.numeric(as.vector(datA2/1000)),main="Posterior prob.: 0.97555",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: A2 (RHGn: Nzime)",cex.names=1)
barplot(as.numeric(as.vector(datB2/1000)),main="Posterior prob.: 0.9703333",ylab="Frequency",xlab="Type of migration",
        names.arg =c("pulse","continuous"),col=col5[1:2],sub="Observed: B2 (RHGn: Ba.Konjo)",cex.names=1)
dev.off()







