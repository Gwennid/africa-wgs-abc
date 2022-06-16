### 20210521
### Gwenna Breton
### Goal: Prepare input files for model choice with ABC-RF. QC, most parameters constant, low and constant migration rates. Based on code from Paul Verdu (PrepINPUT_ABC_Gwenna_Prelim.R) and on 20200312_Prepare_inputs_for_modelchoice.R
### 1000 repeats per model, 2 models per topology.
#Starting with model 1a.

setwd("PLACEHOLDER")

############## BuildInputFiles

nbSim <- 1000
nbModels <- 2

########################################  Sumstats

sumstatsALL <- read.table("3b/ABC_scenario_3b_2models_1000repeats_sumstats_382sumstatswithprop", header=TRUE) #Replace with name of sumstat file

nbSumstat <- length(colnames(sumstatsALL))
nbSumstat

SumStatFULL1 <- matrix(nrow=(nbModels*nbSim+1), ncol=nbSumstat+1)

SumStatFULL1[1,1] <- as.character("x")
SumStatFULL1[1,(2:(nbSumstat+1))] <- as.character(colnames(sumstatsALL))
SumStatFULL1[(2:(nbModels*nbSim+1)),1] <- c(1:(nbModels*nbSim))

SumStatFULL1[(2:(nbModels*nbSim+1)),(2:(nbSumstat+1))] <- as.matrix(sumstatsALL[1:(nbModels*nbSim),])

write.table(SumStatFULL1, file="3b/ABC_scenario_3b_2models_1000repeats.SUMSTATALL", sep="\t", col.names=F, row.names=F, quote=F) #Replace with output name

#################################### Caution! Delete the "x" in the first row (open the file with nano or vim)

#Model Index Table

ModIndexTEMP <- read.table("3b/ABC_scenario_3b_2models_1000repeats_whichmodel", header=F) #Replace with name of file saying which simulation was performed with which model

ModIndex <- matrix(nrow=nbModels*nbSim+1, ncol=2)
ModIndex[1,] <- as.character(c("x","NA"))
ModIndex[(2:(nbModels*nbSim+1)),1] <- c(1:(nbModels*nbSim))
ModIndex[(2:(nbModels*nbSim+1)),2] <- as.character(ModIndexTEMP[,1])


################## IMPORTANT
##########
##########
##########		Double check the file saying which model code (e.g. "1") corresponds to which model definition (e.g. "pulse migration with possibility of high migration rate").
##########
##########
##########

write.table(ModIndex,file="3b/ABC_scenario_3b_2models_1000repeats.MODINDEX", quote=F, sep="\t", col.names=F, row.names=F) #Replace with output name


#################################### Caution! Delete the "NA" in the first row (open the file with nano or vim)
