In R
TODO clean this! Copied from compute_pop_het.R

#####
### Compute observed and expected heterozygosity
#####

## 
# Input files with the count of each configuration are prepared via the commands in compute_pop_het.sh

#write.table(file="allchr_allpop_Het_obs_exp_bialsites.txt",t(c("CHR","POPNAME","H_O_VARSITES","H_O_ALLSITES","H_E_VARSITES","H_E_ALLSITES","H_E_UNBIASED_VARSITES","H_E_UNBIASED_ALLS
ITES","NOVAR_IN_POP","N_TOT_IN_POP")),
#            col.names = FALSE,row.names = FALSE)
write.table(file="allchr_noHGDP_Het_obs_exp_bialsites.txt",t(c("CHR","POPNAME","H_O_VARSITES","H_O_ALLSITES","H_E_VARSITES","H_E_ALLSITES","H_E_UNBIASED_VARSITES","H_E_UNBIASED_ALLSI
TES","NOVAR_IN_POP","N_TOT_IN_POP")),
            col.names = FALSE,row.names = FALSE)
for (CHR in c(1:22)) {
  # for (POPNAME in c("Baka","Nzime","BaKola","Ngumba","AkaMbati","BaKiga","BaTwa","Nsua","BaKonjo","Karretjiepeople","GuiandGana","Juhoansi","Nama","Xun","BantuHerero",
  #                   "BantuKenya","BantuTswana","Biaka","Esan","Gambian","Igbo","Juhoansi_comp","Khomani","Lemande","Luhya","Luo","Mandenka","Maasai","Mbuti",
  #                   "Mende","Mozabite","Saharawi","Yoruba","DaiChinese","Karitiana","Papuan","French","CEU","Coloured","SothoSpeakers","XhosaSpeakers")) {
  for (POPNAME in c("Juhoansi_comp_noHGDP","Dinka_noHGDP","French_noHGDP","Papuan_noHGDP","DaiChinese_noHGDP","Yoruba_noHGDP","Mbuti_noHGDP","Mandenka_noHGDP","Karitiana_noHGDP")) {
    # READ INPUT FILE
    bial_count <- read.table(file=paste("chr",CHR,"_",POPNAME,"_bial_count_genotypes_summary.txt",sep=""),header=TRUE)
    count = read.table(file=paste("chr",CHR,"_",POPNAME,"_count.txt",sep=""))
    # CHECK that HOM_REF + HOM_ALT + HET is always equal to the same thing and set that value to nsam
    if (length(unique(bial_count$HOM_REF + bial_count$HOM_ALT + bial_count$HET)) == 1) {
      nsam = unique(bial_count$HOM_REF + bial_count$HOM_ALT + bial_count$HET)
    } else {
      print(paste("We have an issue! Check the data for chr ",CHR," and population ",POPNAME, sep=""))
    }
    # OBSERVED HET
    ## for each configuration
    H_obs_persite <- bial_count$HET/nsam
    ## on average
    n_tot_obs <- sum(bial_count$X1) #total number of observation
    H_obs <- sum(H_obs_persite*bial_count$X1)/n_tot_obs
    # H_obs can be divided by the total number of invariant sites (in that population) where there is no missingness
    n_novar <- count[2,2] + count[4,2] + count[5,2]
    n_tot <- n_tot_obs+n_novar
    H_obs_allsites <- ((H_obs*n_tot_obs)+0*n_novar)/n_tot
    # EXPECTED HET
    p_ref = (bial_count$HOM_REF + bial_count$HET/2)/nsam
    H_exp_persite = 2*p_ref*(1-p_ref)
    H_exp <- sum(H_exp_persite*bial_count$X1)/n_tot_obs
    H_exp_allsites <- ((H_exp*n_tot_obs)+0*n_novar)/(n_tot_obs+n_novar)
    # UNBIASED EXPECTED HET
    H_exp_unbiased <- (nsam*2)/(nsam*2-1)*H_exp
    H_exp_allsites_unbiased <- (nsam*2)/(nsam*2-1)*H_exp_allsites
    # WRITE OUT
    # write.table(file="allchr_allpop_Het_obs_exp_bialsites.txt",t(c(CHR,POPNAME,H_obs,H_obs_allsites,H_exp,H_exp_allsites,H_exp_unbiased,H_exp_allsites_unbiased,n_novar,n_tot)),append = TRUE,
    #             col.names = FALSE,row.names = FALSE)
    write.table(file="allchr_noHGDP_Het_obs_exp_bialsites.txt",t(c(CHR,POPNAME,H_obs,H_obs_allsites,H_exp,H_exp_allsites,H_exp_unbiased,H_exp_allsites_unbiased,n_novar,n_tot)),append = TRUE,
                col.names = FALSE,row.names = FALSE)
  }
}

#average over the genome
data <- read.table(file="allchr_allpop_Het_obs_exp_bialsites.txt",header=TRUE)
#data <- read.table(file="allchr_noHGDP_Het_obs_exp_bialsites.txt",header=TRUE)
#file with all populations: allpop_avg_het_1-22
write.table(file="allpop_avg_het_1-22",t(c("POPNAME","AVG_EXPECTED_HET","AVG_UNBIASED_EXPECTED_HET","N_ACCESSIBLE_SITES")),row.names = FALSE, col.names = FALSE)
for (POP in c("Baka","Nzime","BaKola","Ngumba","AkaMbati","BaKiga","BaTwa","Nsua","BaKonjo","Karretjiepeople","GuiandGana","Juhoansi","Nama","Xun","BantuHerero",
                    "BantuKenya","BantuTswana","Biaka","Esan","Gambian","Igbo","Juhoansi_comp","Khomani","Lemande","Luhya","Luo","Mandenka","Maasai","Mbuti",
                    "Mende","Mozabite","Saharawi","Yoruba","DaiChinese","Karitiana","Papuan","French","CEU","Coloured","SothoSpeakers","XhosaSpeakers")) {
#for (POP in c("Juhoansi_comp_noHGDP","Dinka_noHGDP","French_noHGDP","Papuan_noHGDP","DaiChinese_noHGDP","Yoruba_noHGDP","Mbuti_noHGDP","Mandenka_noHGDP","Karitiana_noHGDP")) {
pop <- data[data$POPNAME==POP,]
avg_het <- sum(pop$H_E_ALLSITES*pop$N_TOT_IN_POP)/sum(pop$N_TOT_IN_POP)
avg_unbiased_het <- sum(pop$H_E_UNBIASED_ALLSITES*pop$N_TOT_IN_POP)/sum(pop$N_TOT_IN_POP)
write.table(file="allpop_avg_het_1-22",append=TRUE,t(c(POP,avg_het,avg_unbiased_het,sum(pop$N_TOT_IN_POP))),row.names = FALSE, col.names = FALSE)
}

##### X CHROMOSOME

write.table(file="chrX_allpop_Het_exp_bialsites.txt",t(c("CHR","POPNAME","H_E_VARSITES","H_E_ALLSITES","NOVAR_IN_POP","N_TOT_IN_POP")),
#             col.names = FALSE,row.names = FALSE)
# # N_TOT_IN_POP is the number of sites without missing data in the pop which are either non variant or biallelic
# for (POPNAME in c("Baka","Nzime","BaKola","Ngumba","AkaMbati","BaKiga","BaTwa","Nsua","BaKonjo","Karretjiepeople","GuiandGana","Juhoansi","Xun","BantuHerero",
#               "BantuKenya","BantuTswana","Biaka","Esan","Gambian","Igbo","Juhoansi_comp","Khomani","Lemande","Luhya","Luo","Mandenka","Maasai","Mbuti",
#               "Mende","Mozabite","Saharawi","Yoruba","Karitiana","Papuan","French","CEU","Coloured","SothoSpeakers")) {
#   bial_count <- read.table(file=paste("chrX_",POPNAME,"_bial_count_genotypes_summary.txt",sep=""),header=TRUE)
#   count = read.table(file=paste("chrX_",POPNAME,"_count.txt",sep=""))
#   # CHECK that HOM_REF + HOM_ALT + HET is always equal to the same thing and set that value to nsam
#   if (length(unique(2*bial_count$DIPLOID_HOM_REF + 2*bial_count$DIPLOID_HOM_ALT + 2*bial_count$HET + bial_count$HAPLOID_REF + bial_count$HAPLOID_ALT)) == 1) {
#     nsam = unique(2*bial_count$DIPLOID_HOM_REF + 2*bial_count$DIPLOID_HOM_ALT + 2*bial_count$HET + bial_count$HAPLOID_REF + bial_count$HAPLOID_ALT)
#     #caution! Here nsam is the number of chromosome not the number of individuals.
#   } else {
#     print(paste("We have an issue! Check the data for chr X and population ",POPNAME, sep=""))
#   }
#   # OBSERVED HET : not calculated (sample size for females is too low and very variable across populations; and I do not want to make 'fake females').
#   n_tot_obs <- sum(bial_count$X1) #total number of observation TODO CHECK THAT ONE
#   n_novar <- count[2,2] + count[4,2] + count[5,2] # TODO CHECK THAT ONE
#   n_tot <- n_tot_obs+n_novar  # TODO CHECK THAT ONE
#   # EXPECTED HET
#   p_ref = (bial_count$DIPLOID_HOM_REF * 2 + bial_count$HET + bial_count$HAPLOID_REF)/nsam
#   H_exp_persite = 2*p_ref*(1-p_ref)
#   H_exp <- sum(H_exp_persite*bial_count$X1)/n_tot_obs
#   H_exp_allsites <- ((H_exp*n_tot_obs)+0*n_novar)/n_tot
#   # WRITE TO FILE
#   write.table(file="chrX_allpop_Het_exp_bialsites.txt",t(c("X",POPNAME,H_exp,H_exp_allsites,n_novar,n_tot)),
#               col.names = FALSE,row.names = FALSE, append=TRUE)
# }
