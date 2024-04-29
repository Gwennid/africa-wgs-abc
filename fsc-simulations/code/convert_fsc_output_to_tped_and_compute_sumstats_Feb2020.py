#20200227
#Gwenna Breton
#Changes compared to the script used in December2019_secondtry: private alleles fixed (loops over the 100 windows); fix for Tajima's D (remove windows with no variants); do not output the number of non variable sites; work for output of fsc with --dnatosnp 0 option (2s and 3s in biallelic sites are turned to 1s); private alleles: compute average and stdev of the minor allele frequency (not of the derived allele); changed the structure of pos_list to enable computation of e.g. SFS (see fsc_totped_and_python_sumstats_titv0.33 jupyter notebook); added basic implementation of Watterson's theta; added SFS (the proportions); changed the stats related to number of homozygous sites per population to enable easier comparison with observed data; added 'contig' to the name of the windows and sites (first two columns of the TPED); changed to 0 and 1 in the TPED as well.

import sys
import numpy as np
import summary_statistics as ss

filename = sys.argv[1]

fsc_out=open(filename + ".arp")
fsc_out_content = fsc_out.readlines()
CHR = []
ID = []
GEN = []
POS = []
markers_by_chr = [0] # (i+1) - i give the number of polymorphic positions in segment i+1.
count = 0
for i in range(1,101):
    pos = fsc_out_content[19+(i-1)*2].rstrip("\n").strip("#").split(", ")
    count = count + len(pos)
    markers_by_chr.append(count)
    for j in pos:
        CHR.append("contig" + str(i))
        ID.append("contig" + str(i) + ":" + str(j)) 
        GEN.append("0")
        POS.append(str(j))

byhaplotype = [CHR,ID,GEN,POS]
for i in range(1,6):
    sample_start = fsc_out_content.index("\t\tSampleName=\"Sample " + str(i) + "\"\n") #find where the i-th sample starts
    for j in range(1,11):
        genotype = fsc_out_content[sample_start + 2 + j].strip().split("\t")
        byhaplotype.append(list(genotype[2]))
byhaplotypenp = np.array(byhaplotype)

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
#pos_list will be the same for all (because I am not excluding hom_ref or hom_der at this stage).
bi_chr = []
multi_chr = []
with open(filename + ".tped", 'w') as filehandle:
    for i in range(0,100):
        bi = 0
        multi = 0
        CHR = np.zeros([50,(markers_by_chr[i+1]-markers_by_chr[i])],dtype=int)
        BIL = np.zeros([1,(markers_by_chr[i+1]-markers_by_chr[i])],dtype=int)
        notbial = []
        k = 0
        for j in range(markers_by_chr[i],markers_by_chr[i+1]):
            if len(set(byhaplotypenp[:,j][4:])) == 2:
                a = byhaplotypenp[:,j][4:]
		b = [1 if x!='0' else 0 for x in a] #change 2 or 3 to 1
		CHR[:,k] = b
                #CHR[:,k] = [1 if x!='0' else 0 for x in a] #change 2 or 3 to 1
                BIL[:,k] = byhaplotypenp[:,j][3]
                filehandle.writelines("%s " % item for item in byhaplotypenp[:,j][0:4])
		filehandle.writelines("%s " % item for item in b)
                filehandle.write("\n")
                bi = bi + 1
            elif len(set(byhaplotypenp[:,j][4:])) == 1:
                notbial.append(k)
            elif len(set(byhaplotypenp[:,j][4:])) > 2:
                multi = multi + 1
                notbial.append(k)
            k = k+1
        CHR_bial = np.delete(CHR,notbial,axis=1)       
        BIL_bial = np.delete(BIL,notbial,axis=1)
        count = np.sum(CHR_bial, axis=0)
        pos_list.append(BIL_bial[0])
        count_list.append(count)
        hap_list.append(CHR_bial)
        bi_chr.append(bi)
        multi_chr.append(multi)
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
filehandle.close()

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
for i in range(0,100):
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
for i in range(0,100):
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
for i in range(0,100):
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
for i in range(0,100):
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
for i in range(0,100):
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
for i in range(0,100):
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

#Edit 20200227
#Additional summary statistics

##Watterson's theta
#Comment: this is a very basic implementation. It is just the number of biallelic sites divided by the harmonic number for the sample size. TODO think about how to do it properly. Should it be divided by the number of sites at some point? Take stdev across windows? etc.
harmonic_number_9 = 1.0 + 1.0/2.0 + 1.0/3.0 + 1.0/4.0 + 1.0/5.0 + 1.0/6.0 + 1.0/7.0 + 1.0/8.0 + 1.0/9.0
wat_theta_pop1 = float(bi_tot_pop1) / harmonic_number_9
wat_theta_pop2 = float(bi_tot_pop2) / harmonic_number_9 
wat_theta_pop3 = float(bi_tot_pop3) / harmonic_number_9 
wat_theta_pop4 = float(bi_tot_pop4) / harmonic_number_9 
wat_theta_pop5 = float(bi_tot_pop5) / harmonic_number_9 

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
sumstats = open(filename + ".sumstats", 'w')
sumstats.write(str(sum(bi_chr)) + "\t" + str(sum(multi_chr)) + "\t" +
               str(bi_tot_pop1) + "\t" + str(bi_tot_pop2) + "\t" + str(bi_tot_pop3) + "\t" + str(bi_tot_pop4) + "\t" + str(bi_tot_pop5) + "\t" +
               str(hom_pop1) + "\t" + str(hom_pop2) + "\t" + str(hom_pop3) + "\t" + str(hom_pop4) + "\t" + str(hom_pop5) + "\t" +
               str(hom_0_prop_pop1) + "\t" + str(hom_0_prop_pop2) + "\t" + str(hom_0_prop_pop3) + "\t" + str(hom_0_prop_pop4) + "\t" + str(hom_0_prop_pop5) + "\t" +
               str(sumstats_pop1['HET']) + "\t" + str(sumstats_pop2['HET']) + "\t" + str(sumstats_pop3['HET']) + "\t" + str(sumstats_pop4['HET']) + "\t" + str(sumstats_pop5['HET']) + "\t" +
               str(sumstats_pop1['HET_std']) + "\t" + str(sumstats_pop2['HET_std']) + "\t" + str(sumstats_pop3['HET_std']) + "\t" + str(sumstats_pop4['HET_std']) + "\t" + str(sumstats_pop5['HET_std']) + "\t" +
               str(sumstats_pop1['tajD']) + "\t" + str(sumstats_pop2['tajD']) + "\t" + str(sumstats_pop3['tajD']) + "\t" + str(sumstats_pop4['tajD']) + "\t" + str(sumstats_pop5['tajD']) + "\t" +
               str(sumstats_pop1['tajD_std']) + "\t" + str(sumstats_pop2['tajD_std']) + "\t" + str(sumstats_pop3['tajD_std']) + "\t" + str(sumstats_pop4['tajD_std']) + "\t" + str(sumstats_pop5['tajD_std']) + "\t" +
               str(n_private_pop1) + "\t" + str(n_private_pop2) + "\t" + str(n_private_pop3) + "\t" + str(n_private_pop4) + "\t" + str(n_private_pop5) + "\t" +
               str(mean_1) + "\t" + str(mean_2) + "\t" + str(mean_3) + "\t" + str(mean_4) + "\t" +
               str(mean_5) + "\t" + str(var_1) + "\t" + str(var_2) + "\t" + str(var_3) + "\t" + str(var_4) + "\t" + str(var_5) + "\t")
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
sumstats.write(str(wat_theta_pop1) + "\t" + str(wat_theta_pop2) + "\t" + str(wat_theta_pop3)  + "\t" + str(wat_theta_pop4)  + "\t" + str(wat_theta_pop5))
sumstats.close()



