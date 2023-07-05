# 2023-07-05
# Gwenna Breton
# Goal: summarize variant counts information by populations. Should be run on outputs from variant_counts.sh
#Based on /Users/gwennabreton/Documents/Previous_work/PhD_work/P2_RHG_KS/sequence_processing/scripts/mapping_and_GATK_final_scripts/HC_BPresolution/countvariants/20200604_summarize_indcounts_bypop.R

#Preliminary: individual ID were updated to more explicit names in file 25KS.48RHG.104comp.HCBP.1-22.recalSNP99.9.recalINDEL99.0.FAIL1FAIL2FAIL3.reheaded.dbsnp156.some_detail_metrics.AVG using the script africa-wgs-abc/descriptive-analyses/sed_old_new_ind_ID.sh 

##Read in the information
setwd("/Users/gwennabreton/Documents/Previous_work/PhD_work/P2_RHG_KS/writing/ms/africa-wgs-abc/results/variant-counts")
counts <- read.table("25KS.48RHG.104comp.HCBP.1-22.recalSNP99.9.recalINDEL99.0.FAIL1FAIL2FAIL3.reheaded.dbsnp156.some_detail_metrics.AVG_newnames",header=TRUE)

##Add information about population etc
info <- read.table(file="../../descriptive-analyses/25KSP_49PNP_105comp_info.txt",header=TRUE)

counts2 <- counts[1:177,]
info2 <- info[info$ID != "PNP061",]
info3 <- info2[info2$ID != "SGDPJUH2",]

##Ensure the same order
counts2r <- counts2[order(match(counts2[,1],info3[,1])),]

POP <- c("ba.kiga","ba.kola","ba.konjo","ba.twa","baka","bantuh","bantuk","bantut","biaka","biakambati","ceu","coloured","dai",
         "dinka","esan","french","gambian","guigana","igbo","juhoansi","juhoansi_comp","karitiana","karretjie","khomani","kongo","lemande","luhya","luo",
         "maasai","mandenka","mandinka","mbuti","mende","mozabite","nama","ngumba","nsua","nzime","papuan","saharawi","somali","sotho",
         "xhosa","xun","yoruba","zulu")

mat <- matrix(0,length(POP),9)

for (i in (1:length(POP))) {
  pop <- POP[i]
  snp_mean <- mean(counts2r$TOTALSNP[info3$POP==pop])
  snp_sd <- sd(counts2r$TOTALSNP[info3$POP==pop])
  cov_mean <- mean(info3$COVERAGE[info3$POP==pop])
  cov_sd <- sd(info3$COVERAGE[info3$POP==pop])
  novel_mean <- mean(counts2r$TOTALSNP[info3$POP==pop]-counts2r$NUM_IN_DB_SNP[info3$POP==pop])
  novel_sd <- sd(counts2r$TOTALSNP[info3$POP==pop]-counts2r$NUM_IN_DB_SNP[info3$POP==pop])
  indel_mean <- mean(counts2r$TOTALINDEL[info3$POP==pop])
  indel_sd <- sd(counts2r$TOTALINDEL[info3$POP==pop])
  mat[i,] <- c(pop,snp_mean,snp_sd,cov_mean,cov_sd,novel_mean,novel_sd,indel_mean,indel_sd)
}

write.table(file="25KS.48RHG.104comp.allfilters.summarybypop",mat,
            col.names = c("POP","TOTALSNP_mean","TOTALSNP_sd","COV_mean","COV_sd","NOVELSNP_mean","NOVELSNP_sd","TOTALINDEL_mean","TOTALINDEL_sd"),
            row.names=FALSE)
#Comment: the standard deviation for groups with a single individual is NA (logically).

#####QUESTION
##There is more code to get mean coverage by dataset. Do I need that?
#####