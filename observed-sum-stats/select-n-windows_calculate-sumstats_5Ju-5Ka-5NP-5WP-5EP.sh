#!/bin/bash -l
# NOTE the -l flag!
#SBATCH -J TPED-and-observed_sumstats_5Ju-5Ka-5NP-5WP-5EP
#SBATCH -o TPED-and-observed_sumstats_5Ju-5Ka-5NP-5WP-5EP.output
#SBATCH -e TPED-and-observed_sumstats_5Ju-5Ka-5NP-5WP-5EP.output
#SBATCH --mail-user gwenna.breton@gu.se
#SBATCH --mail-type=FAIL,END
#SBATCH -t 24:0:0
#SBATCH -A uppmax2025-2-467
#SBATCH -p core -n 1

### 20251231
### Gwenna Breton
### Goals:
### 1. Select n windows (strategy 2, patchwork windows with modified coordinates).
### 2. Compute observed summary statistics.
### Revisions round 1 for Nature Comm submission, 3rd reviewer point 5.
### Based on 210615_sumstats_strategy2_5Ju-5Na-5Nz-5Baka-5BaT_withRcode.sh and 251119_make-tped-all1.2Mbpwindows_5Ju-5Na-5NP-5WP-5EP.sh (could be any other combination of populations).
### Modified from 251231_select-724-windows_calculate-sumstats_5Ju-5Ka-5NP-5WP-5EP.sh for inclusion in git repo.
### Runs on Uppmax

#
# 0. Set-up
#

module load bioinfo-tools bcftools/1.10
module load bioinfo-tools tabix/0.2.6
module load bioinfo-tools vcftools/0.1.13
module load bioinfo-tools python/2.7.15
module load bioinfo-tools plink/1.90b4.9

pop="5Ju-5Ka-5NP-5WP-5EP" #Order: NK SK NP WP EP
ind="KSP103,KSP105,KSP106,KSP111,KSP116,KSP062,KSP063,KSP065,KSP067,KSP069,PLACEHOLDER" #Order: NK SK WP EP NP
nwindows="100"

folder_windows=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/HC_BPresolution/3maskrecal.realn/allsites/3_geno01_hwefiltering/DIR_acc_map_indel_masked
folder_100=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/observed_sum_stats/Second_trial/strategy2_patchworkwindows/strategy2_patchworkwindows_${nwindows}_1.2Mbpwindows
folder_VCF=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/HC_BPresolution/3maskrecal.realn/allsites/3_geno01_hwefiltering
folder_working=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/observed_sum_stats/Second_trial/strategy2_patchworkwindows/251211_${pop}_${nwindows}_1.2Mbpwindows/
folder_code_obs=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/observed_sum_stats/Second_trial/code
folder_code_simu=/crex/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/Second_trial/code
infile=${folder_working}${pop}.${nwindows}.1.2Mbpwindows.bial.nomiss.renamed.newpos

mkdir ${folder_working}
cd ${folder_working}

#
# 1. Select 100 windows, make TPED
#

while IFS= read -r line 
do
    #Make a VCF with a subset of positions and a subset of individuals
    chr=$(echo $line | cut -f1 -d" ") #something like "chr1"
    int=$(echo $line | cut -f2 -d" ") #something like "1"
    contig=$(echo $line | cut -f6 -d" ") #something like "contig1"
    CHR=$(echo $chr | sed 's/chr//')
    vcf=${folder_VCF}/25KS.48RHG.104comp.HCBP.${CHR}.recalSNP99.9.recalINDEL99.0.FAIL1FAIL2FAIL3.reheaded.PASSSNP.H.vcf.gz
    intervalfile=${folder_windows}/1Mbpwindows/${chr}_${int}
    bcftools view -R ${intervalfile} -s $ind --min-ac 1 -O z ${vcf} -o ${pop}.PASSSNP.${contig}.vcf.gz
    vcf=${folder_working}${pop}.PASSSNP.${contig}
    tabix -f ${vcf}.vcf.gz
    vcftools --gzvcf ${vcf}.vcf.gz --max-missing 1.0 --plink-tped --out ${vcf}.bial.nomiss
    #Rename the variants to chr:pos and rename the chr to contig
    awk '{print $1 ":" $4}' < ${vcf}.bial.nomiss.tped > ${vcf}.bial.nomiss.tped.col2
    cut -f3,4 < ${vcf}.bial.nomiss.tped > ${vcf}.bial.nomiss.tped.col3-4
    cut -f5-54 < ${vcf}.bial.nomiss.tped > ${vcf}.bial.nomiss.tped.col5-54
    nrow=$(wc -l ${vcf}.bial.nomiss.tped | cut -f1 -d" ")
    for i in $(eval echo "{1..${nrow}}"); do echo "${contig}" >> ${vcf}.bial.nomiss.tped.col1 ; done
    paste ${vcf}.bial.nomiss.tped.col1 ${vcf}.bial.nomiss.tped.col2 ${vcf}.bial.nomiss.tped.col3-4 ${vcf}.bial.nomiss.tped.col5-54 > ${vcf}.bial.nomiss.renamed.tped
    rm *col1 *col2 *col3-4
    rm ${vcf}.bial.nomiss.tped
    rm ${vcf}.bial.nomiss.log
    #Change the positions to make them continuous.
    python ${folder_code_obs}/250722_strategy2_real_fake_positions.py ${chr} ${int} ${pop} ${folder_working} ${contig}
    cut -f1-3 ${vcf}.bial.nomiss.renamed.tped > ${vcf}.bial.nomiss.renamed.tped.col1-3
    paste ${vcf}.bial.nomiss.renamed.tped.col1-3 ${contig}_newposcol4 ${vcf}.bial.nomiss.tped.col5-54 > ${vcf}.bial.nomiss.renamed.newpos.tped
    rm ${vcf}.bial.nomiss.renamed.tped.col1-3 ${vcf}.bial.nomiss.tped.col5-54 ${vcf}.bial.nomiss.renamed.tped ${contig}* 
done < ${folder_100}/100windows

##Join all of the tpeds
cat ${pop}.PASSSNP.contig{1..100}.bial.nomiss.renamed.newpos.tped > ${infile}.tped 2>/dev/null

#
# 2. compute summary statistics
#

## ASD
asd --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --out ${infile} --biallelic
Rscript --vanilla ${folder_code_simu}/compute_statistics_based_on_ASD_forrackham.R ${infile} ${folder_working}
#output: ${infile}.asd.dist.sumstats

## Number of bi and multi-allelic variants: can be done directly on the vcf file:
rm ${pop}_multi ${pop}_bi ALT_sorted_counts_${pop}
for contig in {1..100}; do
zgrep -v "#" ${pop}.PASSSNP.contig${contig}.vcf.gz | cut -f5 | sort | uniq -c > ALT_sorted_counts_${pop}
awk 'length($2)>1 {print $1}' ALT_sorted_counts_${pop} >> ${pop}_multi
awk 'length($2)==1 {print $1}' ALT_sorted_counts_${pop} >> ${pop}_bi
done
rm ALT_sorted_counts_${pop}
wc -l ${pop}_bi
##Sanity check: the file ${pop}_bi should have 400 lines (4 different bases * 100 files). TODO check that and/or error message (if a file is not found).

#Write out the total number of biallelic sites and of multiallelic sites
Rscript ${folder_code_obs}/get_bi-multi-tot_n_1.2Mbpwindows.R ${pop}
#Output: ${pop}.100.1.2Mbpwindows.bi.multi.tot.sumstats
#Obs! The Rscript gets the number of contigs from the number of rows in the "_bi" file (divided by 4), so if there is something odd going on, it won't say "100" in the file name.

## Basic sumstats (which do not depend on ancestral state): Python script.
tr '\t' ' ' < ${infile}.tped > tmp
mv tmp ${infile}.tped

for contig in {1..100} ; do
tr '\t' ' ' < ${pop}.PASSSNP.contig${contig}.bial.nomiss.renamed.newpos.tped > tmp
mv tmp ${pop}.PASSSNP.contig${contig}.bial.nomiss.renamed.newpos.tped
done

python ${folder_code_obs}/251229_convertTPED_and_compute_sumstats_nwindows.py ${pop} ${nwindows}
#Output ${infile}.sumstatsH

## ROH sumstats: plink and awk
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --homozyg --homozyg-snp 200 --homozyg-kb 100 --homozyg-density 20 --homozyg-gap 50 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --allow-extra-chr --out ${infile} --missing-genotype 9

### Number of ROH per population
n_ROH_NK=$(grep NK ${infile}.hom | wc -l|cut -f1)
n_ROH_SK=$(grep SK ${infile}.hom | wc -l|cut -f1)
n_ROH_WP=$(grep WP ${infile}.hom | wc -l|cut -f1)
n_ROH_EP=$(grep EP ${infile}.hom | wc -l|cut -f1)
n_ROH_NP=$(grep NP ${infile}.hom | wc -l|cut -f1)

### Number of ROH per class per population
n_ROH_200_NK=$(awk '$1 == "NK" && $9 < 200 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_500_NK=$(awk '$1 == "NK" && $9 > 200 && $9 < 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_500_NK=$(awk '$1 == "NK" && $9 > 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_SK=$(awk '$1 == "SK" && $9 < 200 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_500_SK=$(awk '$1 == "SK" && $9 > 200 && $9 < 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_500_SK=$(awk '$1 == "SK" && $9 > 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_WP=$(awk '$1 == "WP" && $9 < 200 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_500_WP=$(awk '$1 == "WP" && $9 > 200 && $9 < 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_500_WP=$(awk '$1 == "WP" && $9 > 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_EP=$(awk '$1 == "EP" && $9 < 200 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_500_EP=$(awk '$1 == "EP" && $9 > 200 && $9 < 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_500_EP=$(awk '$1 == "EP" && $9 > 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_NP=$(awk '$1 == "NP" && $9 < 200 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_200_500_NP=$(awk '$1 == "NP" && $9 > 200 && $9 < 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)
n_ROH_500_NP=$(awk '$1 == "NP" && $9 > 500 {print $9}' < ${infile}.hom | wc -l|cut -f1)

## Average length of ROH in each class (per population)
avg_L_ROH_200_NK=$(awk '$1 == "NK" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_500_NK=$(awk '$1 == "NK" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_500_NK=$(awk '$1 == "NK" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_SK=$(awk '$1 == "SK" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_500_SK=$(awk '$1 == "SK" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_500_SK=$(awk '$1 == "SK" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_WP=$(awk '$1 == "WP" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_500_WP=$(awk '$1 == "WP" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_500_WP=$(awk '$1 == "WP" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_EP=$(awk '$1 == "EP" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_500_EP=$(awk '$1 == "EP" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_500_EP=$(awk '$1 == "EP" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_NP=$(awk '$1 == "NP" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_200_500_NP=$(awk '$1 == "NP" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)
avg_L_ROH_500_NP=$(awk '$1 == "NP" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < ${infile}.hom)

## Total ROH length per individual : average and variance in the population
avg_L_ROH_NK=$(awk '$1 == "NK" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < ${infile}.hom.indiv)
var_L_ROH_NK=$(awk '$1 == "NK" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < ${infile}.hom.indiv)
avg_L_ROH_SK=$(awk '$1 == "SK" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < ${infile}.hom.indiv)
var_L_ROH_SK=$(awk '$1 == "SK" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < ${infile}.hom.indiv)
avg_L_ROH_WP=$(awk '$1 == "WP" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < ${infile}.hom.indiv)
var_L_ROH_WP=$(awk '$1 == "WP" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < ${infile}.hom.indiv)
avg_L_ROH_EP=$(awk '$1 == "EP" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < ${infile}.hom.indiv)
var_L_ROH_EP=$(awk '$1 == "EP" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < ${infile}.hom.indiv)
avg_L_ROH_NP=$(awk '$1 == "NP" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < ${infile}.hom.indiv)
var_L_ROH_NP=$(awk '$1 == "NP" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < ${infile}.hom.indiv)

echo -e ${n_ROH_NK}\\t${n_ROH_SK}\\t${n_ROH_WP}\\t${n_ROH_EP}\\t${n_ROH_NP}\\t${n_ROH_200_NK}\\t${n_ROH_200_SK}\\t${n_ROH_200_WP}\\t${n_ROH_200_EP}\\t${n_ROH_200_NP}\\t${n_ROH_200_500_NK}\\t${n_ROH_200_500_SK}\\t${n_ROH_200_500_WP}\\t${n_ROH_200_500_EP}\\t${n_ROH_200_500_NP}\\t${n_ROH_500_NK}\\t${n_ROH_500_SK}\\t${n_ROH_500_WP}\\t${n_ROH_500_EP}\\t${n_ROH_500_NP}\\t${avg_L_ROH_200_NK}\\t${avg_L_ROH_200_SK}\\t${avg_L_ROH_200_WP}\\t${avg_L_ROH_200_EP}\\t${avg_L_ROH_200_NP}\\t${avg_L_ROH_200_500_NK}\\t${avg_L_ROH_200_500_SK}\\t${avg_L_ROH_200_500_WP}\\t${avg_L_ROH_200_500_EP}\\t${avg_L_ROH_200_500_NP}\\t${avg_L_ROH_500_NK}\\t${avg_L_ROH_500_SK}\\t${avg_L_ROH_500_WP}\\t${avg_L_ROH_500_EP}\\t${avg_L_ROH_500_NP}\\t${avg_L_ROH_NK}\\t${avg_L_ROH_SK}\\t${avg_L_ROH_WP}\\t${avg_L_ROH_EP}\\t${avg_L_ROH_NP}\\t${var_L_ROH_NK}\\t${var_L_ROH_SK}\\t${var_L_ROH_WP}\\t${var_L_ROH_EP}\\t${var_L_ROH_NP} > ${infile}.ROHsumstats

##FST sumstats: plink and bash
# Comment: this takes a few minutes. I don't think there is anything that needs scaling.
cp ${folder_code_simu}/NK_SK .
cp ${folder_code_simu}/NK_WP .
cp ${folder_code_simu}/NK_EP .
cp ${folder_code_simu}/NK_NP .
cp ${folder_code_simu}/SK_WP .
cp ${folder_code_simu}/SK_EP .
cp ${folder_code_simu}/SK_NP .
cp ${folder_code_simu}/WP_EP .
cp ${folder_code_simu}/WP_NP .
cp ${folder_code_simu}/EP_NP .

plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.NK_SK --fst --within NK_SK
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.NK_WP --fst --within NK_WP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.NK_EP --fst --within NK_EP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.NK_NP --fst --within NK_NP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.SK_WP --fst --within SK_WP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.SK_EP --fst --within SK_EP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.SK_NP --fst --within SK_NP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.WP_EP --fst --within WP_EP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.WP_NP --fst --within WP_NP
plink --tped ${infile}.tped --tfam ${folder_code_simu}/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out ${infile}.EP_NP --fst --within EP_NP

fst_NK_SK=$(grep Mean ${infile}.NK_SK.log | cut -f2 -d":")
fst_NK_WP=$(grep Mean ${infile}.NK_WP.log | cut -f2 -d":")
fst_NK_EP=$(grep Mean ${infile}.NK_EP.log | cut -f2 -d":")
fst_NK_NP=$(grep Mean ${infile}.NK_NP.log | cut -f2 -d":")
fst_SK_WP=$(grep Mean ${infile}.SK_WP.log | cut -f2 -d":")
fst_SK_EP=$(grep Mean ${infile}.SK_EP.log | cut -f2 -d":")
fst_SK_NP=$(grep Mean ${infile}.SK_NP.log | cut -f2 -d":")
fst_WP_EP=$(grep Mean ${infile}.WP_EP.log | cut -f2 -d":")
fst_WP_NP=$(grep Mean ${infile}.WP_NP.log | cut -f2 -d":")
fst_EP_NP=$(grep Mean ${infile}.EP_NP.log | cut -f2 -d":")

echo -e ${fst_NK_SK}\\t${fst_NK_WP}\\t${fst_NK_EP}\\t${fst_NK_NP}\\t${fst_SK_WP}\\t${fst_SK_EP}\\t${fst_SK_NP}\\t${fst_WP_EP}\\t${fst_WP_NP}\\t${fst_EP_NP} > ${infile}.fst.sumstats
#Output: ${infile}.fst.sumstats

## Sumstats that depend on the knowledge of the ancestral state.
### Step1: Get a TPED with ancestral information.
#PS made a script which takes as input the merged TPED and outputs another TPED with only the sites for which ancestral information is known and in which 0 <> ancestral allele. The sites which are masked in the apes are included as "good ancestral" (that could be changed later).
#I modified the file manually so the first argument in the CL is the TPED name.

python ${folder_code_obs}/parse_ancestral_derived_in_tped_regions_nwindows.py ${infile}.tped
#Output: ${infile}_ancestral_state.tped

### Step2: Split the resulting TPED into 100 contigs (easier to match my existing Python scripts).
for i in {1..100}; do grep "contig${i} " ${infile}_ancestral_state.tped > ${infile}_ancestral_state.contig${i}.tped ; done

###Step3: Compute summary statistics: Python script.
python ${folder_code_obs}/251229_ancestral_state_dependant_stats_nwindows.py ${pop} ${nwindows}
# Output: ${infile}.ancestralstate.sumstats

# PUT ALL TOGETHER

## Add a header to each of the files
tr '\n' '\t' < ${folder_code_obs}/header_ancestralstatssumstats_onerowperstat.txt | sed 's/$/\n/g' | cat - ${infile}.ancestralstate.sumstats > ${infile}.ancestralstate.sumstatsH
#106 sumstats
tr '\n' '\t' < ${folder_code_obs}/header_ASDsumstats_onerowperstat.txt | sed 's/$/\n/g' | cat - ${infile}.asd.dist.sumstats > ${infile}.asd.dist.sumstatsH
#180 sumstats
tr '\n' '\t' < ${folder_code_obs}/header_basicsumstats_onerowperstat.txt | sed 's/$/\n/g' | sed 's/^\t//g' | cat - ${infile}.sumstats > ${infile}.sumstatsH
#51 sumstats
tr '\n' '\t' < ${folder_code_obs}/header_bi.multi.tot_sumstats_onerowperstat.txt | sed 's/$/\n/g' | cat - ${pop}.${nwindows}.1.2Mbpwindows.bi.multi.tot.sumstats > ${pop}.${nwindows}.1.2Mbpwindows.bi.multi.tot.sumstatsH
#2 sumstats
tr '\n' '\t' < ${folder_code_obs}/header_FSTsumstats_onerowperstat.txt | sed 's/$/\n/g' | cat - ${infile}.fst.sumstats > ${infile}.fst.sumstatsH
#10 sumstats
tr '\n' '\t' < ${folder_code_obs}/header_ROHsumstats_onerowperstat.txt | sed 's/$/\n/g' | cat - ${infile}.ROHsumstats > ${infile}.ROHsumstatsH
#45 sumstats

## Modify some statistics to proportions.
Rscript ${folder_code_obs}/get_proportion_sumstats_nwindows.R ${pop} ${nwindows}
#Output ${pop}.100.1.2Mbpwindows.prop.sumstatsH

tr ' ' '\t' < ${pop}.${nwindows}.1.2Mbpwindows.prop.sumstatsH > tmp; mv tmp ${pop}.${nwindows}.1.2Mbpwindows.prop.sumstatsH

##Create a vector of summary statistics matching the order in the simulations (caution! See the limitation above).
cut -f1-2 ${pop}.${nwindows}.1.2Mbpwindows.bi.multi.tot.sumstatsH > col1-2
cut -f2-11 ${pop}.${nwindows}.1.2Mbpwindows.prop.sumstatsH > col3-12
cut -f12-16 ${infile}.ancestralstate.sumstatsH > col13-17
cut -f12-31 ${infile}.sumstatsH > col18-37
cut -f12-16 ${pop}.${nwindows}.1.2Mbpwindows.prop.sumstatsH > col38-42
cut -f37-46 ${infile}.sumstatsH > col43-52
cut -f17-106 ${infile}.ancestralstate.sumstatsH > col53-142
cut -f47-51 ${infile}.sumstatsH > col143-147
cut -f1-10 ${infile}.fst.sumstatsH > col148-157
cut -f1-45 ${infile}.ROHsumstatsH > col158-202
cut -f1-180 ${infile}.asd.dist.sumstatsH > col203-382

paste col1-2 col3-12 col13-17 col18-37 col38-42 col43-52 col53-142 col143-147 col148-157 col158-202 col203-382 > 251231.${pop}.${nwindows}.1.2Mbpwindows.vectorof382sumstats
rm col*

############################
##Delete intermediate files.
rm ${pop}.PASSSNP.contig*.bial.nomiss.renamed.newpos.tped
rm NK_SK NK_WP NK_EP NK_NP SK_WP SK_EP SK_NP WP_EP WP_NP EP_NP
rm ${infile}_ancestral_state.contig*.tped
rm *nosex
rm *fst
rm *vcf.gz *vcf.gz.tbi
cp ${pop}.PASSSNP.contig1.bial.nomiss.tfam ${pop}.PASSSNP.contig1.bial.nomiss.tfam_keep
rm ${pop}.PASSSNP.contig*.bial.nomiss.tfam
