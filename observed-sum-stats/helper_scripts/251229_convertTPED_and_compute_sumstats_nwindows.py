import sys
import numpy as np
import summary_statistics as ss

pop = sys.argv[1]
n_windows = int(sys.argv[2])
pos_list=[]
count_list=[]
hap_list=[]
count_list_pop1 = []
count_list_pop2 = []
count_list_pop3 = []
count_list_pop4 = []
count_list_pop5 = []
hap_list_pop1 = []
hap_list_pop2 = []
hap_list_pop3 = []
hap_list_pop4 = []
hap_list_pop5 = []
markers_by_chr=[0]
bi = []
bi_chr = []

for i in range(0,n_windows):
	contig=i+1
	filename= pop + ".PASSSNP.contig" + str(contig) + ".bial.nomiss.renamed.newpos.tped"
        tped=open(filename)
        tped_content=tped.readlines()
        bi_chr.append(len(tped_content))
        CHR_bial = np.zeros([50,bi_chr[i]],dtype=int) #here "int" means integer (it has nothing to do with the interval)
        BIL_bial = np.zeros([1,bi_chr[i]],dtype=int)
        for j in range(0,bi_chr[i]):
            GEN=tped_content[j].strip().split(" ")[4:]
            firstallele = set(GEN).pop()
            CHR_bial[:,j] = [0 if x==firstallele else 1 for x in GEN] #genotypes with 0 and 1 instead of ACTG
            BIL_bial[0,j] = tped_content[j].strip().split(" ")[3]
            count = np.sum(CHR_bial, axis=0)
        pos_list.append(BIL_bial[0])
        hap_list.append(CHR_bial)
        count_list.append(count)
        hap_list_pop1.append(CHR_bial[0:10,:])
        hap_list_pop2.append(CHR_bial[10:20,:])
        hap_list_pop3.append(CHR_bial[20:30,:])
        hap_list_pop4.append(CHR_bial[30:40,:])
        hap_list_pop5.append(CHR_bial[40:50,:])
        count_list_pop1.append(np.sum(CHR_bial[0:10,:],axis=0))
        count_list_pop2.append(np.sum(CHR_bial[10:20,:],axis=0))
        count_list_pop3.append(np.sum(CHR_bial[20:30,:],axis=0))
        count_list_pop4.append(np.sum(CHR_bial[30:40,:],axis=0))
        count_list_pop5.append(np.sum(CHR_bial[40:50,:],axis=0))
	tped.close()

#Subset to biallelic sites in each population
hap_list_pop1_bial = []
pos_list_pop1_bial = []
count_list_pop1_bial = []
hap_list_pop2_bial = []
pos_list_pop2_bial = []
count_list_pop2_bial = []
hap_list_pop3_bial = []
pos_list_pop3_bial = []
count_list_pop3_bial = []
hap_list_pop4_bial = []
pos_list_pop4_bial = []
count_list_pop4_bial = []
hap_list_pop5_bial = []
pos_list_pop5_bial = []
count_list_pop5_bial = []
#Population 1
hom_0_tot_pop1 = 0
hom_1_tot_pop1 = 0
bi_tot_pop1 = 0
for i in range(0,n_windows):
    hom_0 = []
    hom_1 = []
    bi = []
    for j in range(0,bi_chr[i]):
        if count_list_pop1[i][j] == 0:
            hom_0.append(j)
        elif count_list_pop1[i][j] == 10:
            hom_1.append(j)
        else:    
            bi.append(j)
    hap_list_pop1_bial.append(hap_list_pop1[i][:,bi])
    pos_list_pop1_bial.append(pos_list[i][bi])
    count_list_pop1_bial.append(np.sum(hap_list_pop1[i][:,bi],axis=0))
    hom_0_tot_pop1 = hom_0_tot_pop1 + len(hom_0)
    hom_1_tot_pop1 = hom_1_tot_pop1 + len(hom_1)
    bi_tot_pop1 = bi_tot_pop1 + len(bi)
#Population 2
hom_0_tot_pop2 = 0
hom_1_tot_pop2 = 0
bi_tot_pop2 = 0
for i in range(0,n_windows):
    hom_0 = []
    hom_1 = []
    bi = []
    for j in range(0,bi_chr[i]):
        if count_list_pop2[i][j] == 0:
            hom_0.append(j)
        elif count_list_pop2[i][j] == 10:
            hom_1.append(j)
        else:    
            bi.append(j)
    hap_list_pop2_bial.append(hap_list_pop2[i][:,bi])
    pos_list_pop2_bial.append(pos_list[i][bi])
    count_list_pop2_bial.append(np.sum(hap_list_pop2[i][:,bi],axis=0))
    hom_0_tot_pop2 = hom_0_tot_pop2 + len(hom_0)
    hom_1_tot_pop2 = hom_1_tot_pop2 + len(hom_1)
    bi_tot_pop2 = bi_tot_pop2 + len(bi)
#Population 3
hom_0_tot_pop3 = 0
hom_1_tot_pop3 = 0
bi_tot_pop3 = 0
for i in range(0,n_windows):
    hom_0 = []
    hom_1 = []
    bi = []
    for j in range(0,bi_chr[i]):
        if count_list_pop3[i][j] == 0:
            hom_0.append(j)
        elif count_list_pop3[i][j] == 10:
            hom_1.append(j)
        else:    
            bi.append(j)
    hap_list_pop3_bial.append(hap_list_pop3[i][:,bi])
    pos_list_pop3_bial.append(pos_list[i][bi])
    count_list_pop3_bial.append(np.sum(hap_list_pop3[i][:,bi],axis=0))
    hom_0_tot_pop3 = hom_0_tot_pop3 + len(hom_0)
    hom_1_tot_pop3 = hom_1_tot_pop3 + len(hom_1)
    bi_tot_pop3 = bi_tot_pop3 + len(bi)
#Population 4
hom_0_tot_pop4 = 0
hom_1_tot_pop4 = 0
bi_tot_pop4 = 0
for i in range(0,n_windows):
    hom_0 = []
    hom_1 = []
    bi = []
    for j in range(0,bi_chr[i]):
        if count_list_pop4[i][j] == 0:
            hom_0.append(j)
        elif count_list_pop4[i][j] == 10:
            hom_1.append(j)
        else:    
            bi.append(j)
    hap_list_pop4_bial.append(hap_list_pop4[i][:,bi])
    pos_list_pop4_bial.append(pos_list[i][bi])
    count_list_pop4_bial.append(np.sum(hap_list_pop4[i][:,bi],axis=0))
    hom_0_tot_pop4 = hom_0_tot_pop4 + len(hom_0)
    hom_1_tot_pop4 = hom_1_tot_pop4 + len(hom_1)
    bi_tot_pop4 = bi_tot_pop4 + len(bi)
#Population 5
hom_0_tot_pop5 = 0
hom_1_tot_pop5 = 0
bi_tot_pop5 = 0
for i in range(0,n_windows):
    hom_0 = []
    hom_1 = []
    bi = []
    for j in range(0,bi_chr[i]):
        if count_list_pop5[i][j] == 0:
            hom_0.append(j)
        elif count_list_pop5[i][j] == 10:
            hom_1.append(j)
        else:    
            bi.append(j)
    hap_list_pop5_bial.append(hap_list_pop5[i][:,bi])
    pos_list_pop5_bial.append(pos_list[i][bi])
    count_list_pop5_bial.append(np.sum(hap_list_pop5[i][:,bi],axis=0))
    hom_0_tot_pop5 = hom_0_tot_pop5 + len(hom_0)
    hom_1_tot_pop5 = hom_1_tot_pop5 + len(hom_1)
    bi_tot_pop5 = bi_tot_pop5 + len(bi)

#Private diversity
#Edit 20200227 - count of minor allele
freq_private_minor_pop1 = []
freq_private_minor_pop2 = []
freq_private_minor_pop3 = []
freq_private_minor_pop4 = []
freq_private_minor_pop5 = []
for i in range(0,n_windows):
    for j in range(0,bi_chr[i]):
        if (count_list_pop1[i][j] != 0) and (count_list_pop2[i][j] == 0) and (count_list_pop3[i][j] == 0) and (count_list_pop4[i][j] == 0) and (count_list_pop5[i][j] == 0):
            if (count_list_pop1[i][j] > 5):
                freq_private_minor_pop1.append(10-count_list_pop1[i][j])
            else:
                freq_private_minor_pop1.append(count_list_pop1[i][j])
        elif (count_list_pop1[i][j] == 0) and (count_list_pop2[i][j] != 0) and (count_list_pop3[i][j] == 0) and (count_list_pop4[i][j] == 0) and (count_list_pop5[i][j] == 0):
            if (count_list_pop2[i][j] > 5):
                freq_private_minor_pop2.append(10-count_list_pop2[i][j])
            else:
                freq_private_minor_pop2.append(count_list_pop2[i][j])
        elif (count_list_pop1[i][j] == 0) and (count_list_pop2[i][j] == 0) and (count_list_pop3[i][j] != 0) and (count_list_pop4[i][j] == 0) and (count_list_pop5[i][j] == 0):
            if (count_list_pop3[i][j] > 5):
                freq_private_minor_pop3.append(10-count_list_pop3[i][j])
            else:
                freq_private_minor_pop3.append(count_list_pop3[i][j])
        elif (count_list_pop1[i][j] == 0) and (count_list_pop2[i][j] == 0) and (count_list_pop3[i][j] == 0) and (count_list_pop4[i][j] != 0) and (count_list_pop5[i][j] == 0):
            if (count_list_pop4[i][j] > 5):
                freq_private_minor_pop4.append(10-count_list_pop4[i][j])
            else:
                freq_private_minor_pop4.append(count_list_pop4[i][j])
        elif (count_list_pop1[i][j] == 0) and (count_list_pop2[i][j] == 0) and (count_list_pop3[i][j] == 0) and (count_list_pop4[i][j] == 0) and (count_list_pop5[i][j] != 0):
            if (count_list_pop5[i][j] > 5):
                freq_private_minor_pop5.append(10-count_list_pop5[i][j])
            else:
                freq_private_minor_pop5.append(count_list_pop5[i][j])

#Summary statistics
#Edit 20200130: adding a step to remove empty arrays.
count_list_pop1_bial_noempty = [x for x in count_list_pop1_bial if x.any()]
res_classic = ss.classical_stats(10, count_list_pop1_bial_noempty)
sumstats_pop1 = dict(zip(['HET', 'HET_std', 'PI', 'PI_std', 'tajD', 'tajD_std'] , res_classic))

count_list_pop2_bial_noempty = [x for x in count_list_pop2_bial if x.any()]
res_classic = ss.classical_stats(10, count_list_pop2_bial_noempty)
sumstats_pop2 = dict(zip(['HET', 'HET_std', 'PI', 'PI_std', 'tajD', 'tajD_std'] , res_classic))

count_list_pop3_bial_noempty = [x for x in count_list_pop3_bial if x.any()]
res_classic = ss.classical_stats(10, count_list_pop3_bial_noempty)
sumstats_pop3 = dict(zip(['HET', 'HET_std', 'PI', 'PI_std', 'tajD', 'tajD_std'] , res_classic))

count_list_pop4_bial_noempty = [x for x in count_list_pop4_bial if x.any()]
res_classic = ss.classical_stats(10, count_list_pop4_bial_noempty)
sumstats_pop4 = dict(zip(['HET', 'HET_std', 'PI', 'PI_std', 'tajD', 'tajD_std'] , res_classic))

count_list_pop5_bial_noempty = [x for x in count_list_pop5_bial if x.any()]
res_classic = ss.classical_stats(10, count_list_pop5_bial_noempty)
sumstats_pop5 = dict(zip(['HET', 'HET_std', 'PI', 'PI_std', 'tajD', 'tajD_std'] , res_classic))

#Private diversity summary statistics
n_private_pop1 = len(freq_private_minor_pop1) #number of sites private to population 1
n_private_pop2 = len(freq_private_minor_pop2)
n_private_pop3 = len(freq_private_minor_pop3)
n_private_pop4 = len(freq_private_minor_pop4)
n_private_pop5 = len(freq_private_minor_pop5)
mean_1 = np.mean(np.array(freq_private_minor_pop1) / 10.0) #average allelic frequency of variants private to population 1
var_1 = np.var(np.array(freq_private_minor_pop1) / 10.0) #variance of the allelic frequency of variants private to population 1
mean_2 = np.mean(np.array(freq_private_minor_pop2) / 10.0)
var_2 = np.var(np.array(freq_private_minor_pop2) / 10.0)
mean_3 = np.mean(np.array(freq_private_minor_pop3) / 10.0)
var_3 = np.var(np.array(freq_private_minor_pop3) / 10.0)
mean_4 = np.mean(np.array(freq_private_minor_pop4) / 10.0)
var_4 = np.var(np.array(freq_private_minor_pop4) / 10.0)
mean_5 = np.mean(np.array(freq_private_minor_pop5) / 10.0)
var_5 = np.var(np.array(freq_private_minor_pop5) / 10.0)

#Edit 20200309
#Additional summary statistics (the ones which do not require knowledge of ancestral state)

##Watterson's theta
#Comment: this is a very basic implementation. It is just the number of biallelic sites divided by the harmonic number for the sample size. TODO think about how to do it properly. Should it be divided by the number of sites at some point? Take stdev across windows? etc.
harmonic_number_9 = 1.0 + 1.0/2.0 + 1.0/3.0 + 1.0/4.0 + 1.0/5.0 + 1.0/6.0 + 1.0/7.0 + 1.0/8.0 + 1.0/9.0
wat_theta_pop1 = float(bi_tot_pop1) / harmonic_number_9
wat_theta_pop2 = float(bi_tot_pop2) / harmonic_number_9 
wat_theta_pop3 = float(bi_tot_pop3) / harmonic_number_9 
wat_theta_pop4 = float(bi_tot_pop4) / harmonic_number_9 
wat_theta_pop5 = float(bi_tot_pop5) / harmonic_number_9 

##Homozygous sites per population
#Sum of homozygous sites
hom_pop1 = hom_0_tot_pop1 + hom_1_tot_pop1
hom_pop2 = hom_0_tot_pop2 + hom_1_tot_pop2
hom_pop3 = hom_0_tot_pop3 + hom_1_tot_pop3
hom_pop4 = hom_0_tot_pop4 + hom_1_tot_pop4
hom_pop5 = hom_0_tot_pop5 + hom_1_tot_pop5

#Write out everything - comment: less statistics than in the script for the simulations (other statistics will be computed in other scripts).
sumstats = open(pop + "." + str(n_windows) + ".1.2Mbpwindows.bial.nomiss.renamed.newpos.sumstats", 'w')
sumstats.write(str(sum(bi_chr)) + "\t" +
               str(bi_tot_pop1) + "\t" + str(bi_tot_pop2) + "\t" + str(bi_tot_pop3) + "\t" + str(bi_tot_pop4) + "\t" + str(bi_tot_pop5) + "\t" +
               str(hom_pop1) + "\t" + str(hom_pop2) + "\t" + str(hom_pop3) + "\t" + str(hom_pop4) + "\t" + str(hom_pop5) + "\t" +
               str(sumstats_pop1['HET']) + "\t" + str(sumstats_pop2['HET']) + "\t" + str(sumstats_pop3['HET']) + "\t" + str(sumstats_pop4['HET']) + "\t" + str(sumstats_pop5['HET']) + "\t" +
               str(sumstats_pop1['HET_std']) + "\t" + str(sumstats_pop2['HET_std']) + "\t" + str(sumstats_pop3['HET_std']) + "\t" + str(sumstats_pop4['HET_std']) + "\t" + str(sumstats_pop5['HET_std']) + "\t" +
               str(sumstats_pop1['tajD']) + "\t" + str(sumstats_pop2['tajD']) + "\t" + str(sumstats_pop3['tajD']) + "\t" + str(sumstats_pop4['tajD']) + "\t" + str(sumstats_pop5['tajD']) + "\t" +
               str(sumstats_pop1['tajD_std']) + "\t" + str(sumstats_pop2['tajD_std']) + "\t" + str(sumstats_pop3['tajD_std']) + "\t" + str(sumstats_pop4['tajD_std']) + "\t" + str(sumstats_pop5['tajD_std']) + "\t" +
               str(n_private_pop1) + "\t" + str(n_private_pop2) + "\t" + str(n_private_pop3) + "\t" + str(n_private_pop4) + "\t" + str(n_private_pop5) + "\t" +
               str(mean_1) + "\t" + str(mean_2) + "\t" + str(mean_3) + "\t" + str(mean_4) + "\t" +
               str(mean_5) + "\t" + str(var_1) + "\t" + str(var_2) + "\t" + str(var_3) + "\t" + str(var_4) + "\t" + str(var_5) + "\t" + str(wat_theta_pop1) + "\t" + str(wat_theta_pop2) + "\t" + str(wat_theta_pop3)  + "\t" + str(wat_theta_pop4)  + "\t" + str(wat_theta_pop5))
sumstats.close()



