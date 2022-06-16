### 20210521
### Gwenna Breton
### Goal: Prepare input files for model choice with ABC-RF. QC, most parameters constant, low and constant migration rates. Based on code from Paul Verdu (PrepINPUT_ABC_Gwenna_Prelim.R) and on 20200312_Prepare_inputs_for_modelchoice.R
### 500 repeats per model, 2 models per topology.
#Starting with model 1a.

setwd("/home/gwennabreton/Documents/PhD/Projects/Project2_Pygmies_KhoeSan_project/sequence_processing/scripts/ABC/modelchoice_pulseVScontinuous/QC/Constant_parameters/")

############## BuildInputFiles

nbSim <- 1000
nbModels <- 2

########################################  Sumstats

sumstatsALL <- read.table("3b/ABC_scenario_3b_2models_1000repeats_sumstats_382sumstatswithprop", header=TRUE)

nbSumstat <- length(colnames(sumstatsALL))
nbSumstat

SumStatFULL1 <- matrix(nrow=(nbModels*nbSim+1), ncol=nbSumstat+1)

SumStatFULL1[1,1] <- as.character("x")
SumStatFULL1[1,(2:(nbSumstat+1))] <- as.character(colnames(sumstatsALL))
SumStatFULL1[(2:(nbModels*nbSim+1)),1] <- c(1:(nbModels*nbSim))

SumStatFULL1[(2:(nbModels*nbSim+1)),(2:(nbSumstat+1))] <- as.matrix(sumstatsALL[1:(nbModels*nbSim),])

write.table(SumStatFULL1, file="3b/ABC_scenario_3b_2models_1000repeats.SUMSTATALL", sep="\t", col.names=F, row.names=F, quote=F)

#################################### ACHTUNG !!! Efface le x (directement dans le fichier, avec nano).

#Model Index Table

ModIndexTEMP <- read.table("3b/ABC_scenario_3b_2models_1000repeats_whichmodel", header=F)

ModIndex <- matrix(nrow=nbModels*nbSim+1, ncol=2)
ModIndex[1,] <- as.character(c("x","NA"))
ModIndex[(2:(nbModels*nbSim+1)),1] <- c(1:(nbModels*nbSim))
ModIndex[(2:(nbModels*nbSim+1)),2] <- as.character(ModIndexTEMP[,1])


################## IMPORTANT
##########
##########
##########		Check file "correspondence_modelnumber.txt" for Model Index vs Def originale des modèles. First three rows: pulse, last three rows: continuous.
##########
##########
##########

write.table(ModIndex,file="3b/ABC_scenario_3b_2models_1000repeats.MODINDEX", quote=F, sep="\t", col.names=F, row.names=F)


#################################### ACHTUNG !!! Efface le NA



