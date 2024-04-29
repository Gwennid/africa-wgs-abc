
module load bioinfo-tools python/2.7.15
module load bioinfo-tools R/3.6.1
module load bioinfo-tools plink/1.90b4.9

folder=YOURFOLDER #TODO update
modelfolder=${folder}/models_and_parameters_sets

#Copy everything to scratch.
cd $SNIC_TMP
mkdir $SNIC_TMP/ASDmatrices
cp ${modelfolder}/${scenario}.def ${modelfolder}/${scenario}.def_header ${modelfolder}/${scenario}_raw.tpl ${folder}/code/convert_fsc_output_to_tped_and_compute_sumstats_Feb2020.py ${folder}/code/summary_statistics.py ${folder}/code/NK_SK_WP_EP_NP.tfam ${folder}/code/NK_SK ${folder}/code/NK_WP ${folder}/code/NK_EP ${folder}/code/NK_NP ${folder}/code/SK_WP ${folder}/code/SK_EP ${folder}/code/SK_NP ${folder}/code/WP_EP ${folder}/code/WP_NP ${folder}/code/EP_NP ${folder}/code/go_check_this_simulation.txt ${folder}/code/file_is_empty.txt ${folder}/code/file_does_not_exist.txt .

for i in $(eval echo "{$start..$end}"); do

#Step 0: select a vector of parameter values and create a temporary folder for that simulation.
mkdir $SNIC_TMP/tmp_${scenario}_${i}
cd $SNIC_TMP/tmp_${scenario}_${i}/
k=$(($i+1))
head -$k $SNIC_TMP/${scenario}.def | tail -n1 > $SNIC_TMP/tmp_${scenario}_${i}/param_vector
cat $SNIC_TMP/${scenario}.def_header $SNIC_TMP/tmp_${scenario}_${i}/param_vector > $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_$i.def
cp $SNIC_TMP/${scenario}_raw.tpl $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}.tpl

#Step 1: run fsc
fsc26 -t $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}.tpl -n 1 -f $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}.def --quiet --dnatosnp 0 #the dnatosnp option is to output ancestral state as 0.
seed=$(cut -d":" -f2 seed.txt)

#Step 2: run python script which converts fsc output to tped and to input for Flora Jays scripts and outputs a number of summary statistics
python $SNIC_TMP/convert_fsc_output_to_tped_and_compute_sumstats_Feb2020.py $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1

#Step 3: run asd
asd --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1 --biallelic --missing-tped -9 #the --missing-tped option is to avoid asd considering every 0 as missing!

#Step 4: ROH analysis (plink to create them, awk/bash to extract information)
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --homozyg --homozyg-snp 200 --homozyg-kb 100 --homozyg-density 20 --homozyg-gap 50 --homozyg-window-snp 50 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.05 --allow-extra-chr --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1 --missing-genotype 9

## Number of ROH per population
n_ROH_NK=$(grep NK $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_SK=$(grep SK $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_WP=$(grep WP $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_EP=$(grep EP $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_NP=$(grep NP $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)

## Number of ROH per class per population
n_ROH_200_NK=$(awk '$1 == "NK" && $9 < 200 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_500_NK=$(awk '$1 == "NK" && $9 > 200 && $9 < 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_500_NK=$(awk '$1 == "NK" && $9 > 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_SK=$(awk '$1 == "SK" && $9 < 200 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_500_SK=$(awk '$1 == "SK" && $9 > 200 && $9 < 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_500_SK=$(awk '$1 == "SK" && $9 > 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_WP=$(awk '$1 == "WP" && $9 < 200 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_500_WP=$(awk '$1 == "WP" && $9 > 200 && $9 < 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_500_WP=$(awk '$1 == "WP" && $9 > 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_EP=$(awk '$1 == "EP" && $9 < 200 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_500_EP=$(awk '$1 == "EP" && $9 > 200 && $9 < 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_500_EP=$(awk '$1 == "EP" && $9 > 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_NP=$(awk '$1 == "NP" && $9 < 200 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_200_500_NP=$(awk '$1 == "NP" && $9 > 200 && $9 < 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)
n_ROH_500_NP=$(awk '$1 == "NP" && $9 > 500 {print $9}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom | wc -l|cut -f1)

## Average length of ROH in each class (per population)
avg_L_ROH_200_NK=$(awk '$1 == "NK" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_500_NK=$(awk '$1 == "NK" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_500_NK=$(awk '$1 == "NK" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_SK=$(awk '$1 == "SK" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_500_SK=$(awk '$1 == "SK" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_500_SK=$(awk '$1 == "SK" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_WP=$(awk '$1 == "WP" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_500_WP=$(awk '$1 == "WP" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_500_WP=$(awk '$1 == "WP" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_EP=$(awk '$1 == "EP" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_500_EP=$(awk '$1 == "EP" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_500_EP=$(awk '$1 == "EP" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_NP=$(awk '$1 == "NP" && $9 < 200 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_200_500_NP=$(awk '$1 == "NP" && $9 > 200 && $9 < 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)
avg_L_ROH_500_NP=$(awk '$1 == "NP" && $9 > 500 {total += $9; count ++} END {if (count > 0) {print total/count} else {print 0}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom)

## Total ROH length per individual : average and variance in the population
avg_L_ROH_NK=$(awk '$1 == "NK" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
var_L_ROH_NK=$(awk '$1 == "NK" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
avg_L_ROH_SK=$(awk '$1 == "SK" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
var_L_ROH_SK=$(awk '$1 == "SK" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
avg_L_ROH_WP=$(awk '$1 == "WP" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
var_L_ROH_WP=$(awk '$1 == "WP" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
avg_L_ROH_EP=$(awk '$1 == "EP" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
var_L_ROH_EP=$(awk '$1 == "EP" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
avg_L_ROH_NP=$(awk '$1 == "NP" {count[$1]++; sum[$1]+=$5} END {for(i in count) {m = sum[i]/count[i]; print m}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)
var_L_ROH_NP=$(awk '$1 == "NP" {count[$1]++; sum[$1]+=$5; sumsq[$1]+=$5*$5} END {for(i in count) {m = sum[i]/count[i]; print (sumsq[i] - count[i]*m**2)/(count[i]-1)}}' < $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.hom.indiv)

echo -e ${n_ROH_NK}\\t${n_ROH_SK}\\t${n_ROH_WP}\\t${n_ROH_EP}\\t${n_ROH_NP}\\t${n_ROH_200_NK}\\t${n_ROH_200_SK}\\t${n_ROH_200_WP}\\t${n_ROH_200_EP}\\t${n_ROH_200_NP}\\t${n_ROH_200_500_NK}\\t${n_ROH_200_500_SK}\\t${n_ROH_200_500_WP}\\t${n_ROH_200_500_EP}\\t${n_ROH_200_500_NP}\\t${n_ROH_500_NK}\\t${n_ROH_500_SK}\\t${n_ROH_500_WP}\\t${n_ROH_500_EP}\\t${n_ROH_500_NP}\\t${avg_L_ROH_200_NK}\\t${avg_L_ROH_200_SK}\\t${avg_L_ROH_200_WP}\\t${avg_L_ROH_200_EP}\\t${avg_L_ROH_200_NP}\\t${avg_L_ROH_200_500_NK}\\t${avg_L_ROH_200_500_SK}\\t${avg_L_ROH_200_500_WP}\\t${avg_L_ROH_200_500_EP}\\t${avg_L_ROH_200_500_NP}\\t${avg_L_ROH_500_NK}\\t${avg_L_ROH_500_SK}\\t${avg_L_ROH_500_WP}\\t${avg_L_ROH_500_EP}\\t${avg_L_ROH_500_NP}\\t${avg_L_ROH_NK}\\t${avg_L_ROH_SK}\\t${avg_L_ROH_WP}\\t${avg_L_ROH_EP}\\t${avg_L_ROH_NP}\\t${var_L_ROH_NK}\\t${var_L_ROH_SK}\\t${var_L_ROH_WP}\\t${var_L_ROH_EP}\\t${var_L_ROH_NP} > $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.ROHsumstats

#Step 5: FST analysis (plink to calculate them, bash to extract information?)
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_SK --fst --within $SNIC_TMP/NK_SK
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_WP --fst --within $SNIC_TMP/NK_WP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_EP --fst --within $SNIC_TMP/NK_EP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_NP --fst --within $SNIC_TMP/NK_NP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.SK_WP --fst --within $SNIC_TMP/SK_WP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.SK_EP --fst --within $SNIC_TMP/SK_EP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.SK_NP --fst --within $SNIC_TMP/SK_NP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.WP_EP --fst --within $SNIC_TMP/WP_EP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.WP_NP --fst --within $SNIC_TMP/WP_NP
plink --tped $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.tped --tfam $SNIC_TMP/NK_SK_WP_EP_NP.tfam --allow-extra-chr --missing-genotype 9 --out $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.EP_NP --fst --within $SNIC_TMP/EP_NP

fst_NK_SK=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_SK.log | cut -f2 -d":")
fst_NK_WP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_WP.log | cut -f2 -d":")
fst_NK_EP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_EP.log | cut -f2 -d":")
fst_NK_NP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.NK_NP.log | cut -f2 -d":")
fst_SK_WP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.SK_WP.log | cut -f2 -d":")
fst_SK_EP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.SK_EP.log | cut -f2 -d":")
fst_SK_NP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.SK_NP.log | cut -f2 -d":")
fst_WP_EP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.WP_EP.log | cut -f2 -d":")
fst_WP_NP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.WP_NP.log | cut -f2 -d":")
fst_EP_NP=$(grep Mean $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.EP_NP.log | cut -f2 -d":")
echo -e ${fst_NK_SK}\\t${fst_NK_WP}\\t${fst_NK_EP}\\t${fst_NK_NP}\\t${fst_SK_WP}\\t${fst_SK_EP}\\t${fst_SK_NP}\\t${fst_WP_EP}\\t${fst_WP_NP}\\t${fst_EP_NP} > $SNIC_TMP/tmp_${scenario}_${i}/fst

#Step 6: Write the different items to files. I decided to create one file where everything is together, but also files with only subsets (python sumstats / fst mostly) in case something goes wrong - I think it will make it easier to understand what went wrong.
#The paste commands won't work if the files do not exist. They will work if the files are empty though.
echo -e ${i}\\t${seed} > $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed
if [ -f $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.sumstats ] && [ -f $SNIC_TMP/tmp_${scenario}_${i}/fst ] && [ -f $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.ROHsumstats ] ; then
paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/tmp_${scenario}_${i}/param_vector $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.sumstats $SNIC_TMP/tmp_${scenario}_${i}/fst $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.ROHsumstats >> $SNIC_TMP/${scenario}_paramsets_sumstats_${start}-${end};
else
paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/tmp_${scenario}_${i}/param_vector $SNIC_TMP/go_check_this_simulation.txt >> $SNIC_TMP/${scenario}_paramsets_sumstats_${start}-${end} ;
fi

#For the single sets of sumstats, I will check that the file is not empty
if [ -f $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.sumstats ]; then
	if [ -s $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.sumstats ]; then
		paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.sumstats >> $SNIC_TMP/${scenario}_pythonsumstats_${start}-${end}
	else
		paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/file_is_empty.txt >> $SNIC_TMP/${scenario}_pythonsumstats_${start}-${end}
	fi
else
	paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/file_does_not_exist.txt >> $SNIC_TMP/${scenario}_pythonsumstats_${start}-${end}
fi		

if [ -f $SNIC_TMP/tmp_${scenario}_${i}/fst ]; then
	if [ -s $SNIC_TMP/tmp_${scenario}_${i}/fst ]; then
		paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/tmp_${scenario}_${i}/fst >> $SNIC_TMP/${scenario}_fst_${start}-${end}
	else
		paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/file_is_empty.txt >> $SNIC_TMP/${scenario}_fst_${start}-${end}
	fi
else
	paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/file_does_not_exist.txt >> $SNIC_TMP/${scenario}_fst_${start}-${end}
fi

if [ -f $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.ROHsumstats ]; then
	if [ -s $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.ROHsumstats ]; then
		paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.ROHsumstats >> $SNIC_TMP/${scenario}_ROHsumstats_${start}-${end}
	else
		paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/file_is_empty.txt >> $SNIC_TMP/${scenario}_ROHsumstats_${start}-${end}
	fi
else
	paste $SNIC_TMP/tmp_${scenario}_${i}/iterationandseed $SNIC_TMP/file_does_not_exist.txt >> $SNIC_TMP/${scenario}_ROHsumstats_${start}-${end}
fi

#Step 6: Clean
cp $SNIC_TMP/tmp_${scenario}_${i}/${scenario}_${i}/${scenario}_${i}_1_1.asd.dist $SNIC_TMP/ASDmatrices
cd $SNIC_TMP
rm -r $SNIC_TMP/tmp_${scenario}_${i}/

done

#Once all iterations have run: move everything back.
cp $SNIC_TMP/${scenario}_paramsets_sumstats_${start}-${end} ${folder}/output/sumstats/
cp $SNIC_TMP/${scenario}_pythonsumstats_${start}-${end} ${folder}/output/sumstats/
cp $SNIC_TMP/${scenario}_fst_${start}-${end} ${folder}/output/sumstats/
cp $SNIC_TMP/${scenario}_ROHsumstats_${start}-${end} ${folder}/output/sumstats/
cp $SNIC_TMP/ASDmatrices/* ${folder}/output/ASD/






