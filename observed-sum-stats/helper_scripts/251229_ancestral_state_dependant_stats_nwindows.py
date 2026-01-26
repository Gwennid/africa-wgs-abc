#20200309
#Gwenna Breton
#Goal: take the TPED where 0 and 1 are informative in term of ancestral state (obtained after running PS's script and splitting by contig), format the information so it can be used by FJ's Python scripts (and mine), output summary statistics which require knowledge of the ancestral / derived state (plus some stats to check that everything looks fine). 

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
	k = i+1
        filename= pop + "." + n_windows + ".1.2Mbpwindows.bial.nomiss.renamed.newpos_ancestral_state.contig" + str(k) + ".tped"
        tped=open(filename)
        tped_content=tped.readlines()
        bi_chr.append(len(tped_content))
        CHR_bial = np.zeros([50,bi_chr[i]],dtype=int) #here "int" means integer (it has nothing to do with the interval)
        BIL_bial = np.zeros([1,bi_chr[i]],dtype=int)
        for j in range(0,bi_chr[i]):
		b = tped_content[j].strip().split(" ")[4:]
		CHR_bial[:,j] = b
		BIL_bial[:,k] = tped_content[j].strip().split(" ")[3]
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

##SFS (unfolded)
n_haplo=10
# Allele frequency spectrum
res_afs_pop1 = ss.spatial_histo_fast(pos_list_pop1_bial, count_list_pop1_bial, n_haplo-1)
res_afs_pop2 = ss.spatial_histo_fast(pos_list_pop2_bial, count_list_pop2_bial, n_haplo-1)
res_afs_pop3 = ss.spatial_histo_fast(pos_list_pop3_bial, count_list_pop3_bial, n_haplo-1)
res_afs_pop4 = ss.spatial_histo_fast(pos_list_pop4_bial, count_list_pop4_bial, n_haplo-1)
res_afs_pop5 = ss.spatial_histo_fast(pos_list_pop5_bial, count_list_pop5_bial, n_haplo-1)

##Homozygous sites per population
#Sum of homozygous sites
hom_pop1 = hom_0_tot_pop1 + hom_1_tot_pop1
hom_pop2 = hom_0_tot_pop2 + hom_1_tot_pop2
hom_pop3 = hom_0_tot_pop3 + hom_1_tot_pop3
hom_pop4 = hom_0_tot_pop4 + hom_1_tot_pop4
hom_pop5 = hom_0_tot_pop5 + hom_1_tot_pop5
#Proportion of homozygous sites for the ancestral allele
hom_0_prop_pop1 = float(hom_0_tot_pop1) / float(hom_pop1)
hom_0_prop_pop2 = float(hom_0_tot_pop2) / float(hom_pop2)
hom_0_prop_pop3 = float(hom_0_tot_pop3) / float(hom_pop3)
hom_0_prop_pop4 = float(hom_0_tot_pop4) / float(hom_pop4)
hom_0_prop_pop5 = float(hom_0_tot_pop5) / float(hom_pop5)

#write out everything
sumstats = open(pop + "." + str(n_windows) + ".1.2Mbpwindows.bial.nomiss.renamed.newpos.ancestralstate.sumstats", 'w')
sumstats.write(str(sum(bi_chr)) + "\t" +
               str(bi_tot_pop1) + "\t" + str(bi_tot_pop2) + "\t" + str(bi_tot_pop3) + "\t" + str(bi_tot_pop4) + "\t" + str(bi_tot_pop5) + "\t" +
               str(hom_pop1) + "\t" + str(hom_pop2) + "\t" + str(hom_pop3) + "\t" + str(hom_pop4) + "\t" + str(hom_pop5) + "\t" +
               str(hom_0_prop_pop1) + "\t" + str(hom_0_prop_pop2) + "\t" + str(hom_0_prop_pop3) + "\t" + str(hom_0_prop_pop4) + "\t" + str(hom_0_prop_pop5) + "\t")
sumstats.writelines("%s\t" % item for item in res_afs_pop1[0])
sumstats.writelines("%s\t" % item for item in res_afs_pop2[0])
sumstats.writelines("%s\t" % item for item in res_afs_pop3[0])
sumstats.writelines("%s\t" % item for item in res_afs_pop4[0])
sumstats.writelines("%s\t" % item for item in res_afs_pop5[0])
sumstats.writelines("%s\t" % item for item in res_afs_pop1[1])
sumstats.writelines("%s\t" % item for item in res_afs_pop2[1])
sumstats.writelines("%s\t" % item for item in res_afs_pop3[1])
sumstats.writelines("%s\t" % item for item in res_afs_pop4[1])
sumstats.writelines("%s\t" % item for item in res_afs_pop5[1])
sumstats.close()

