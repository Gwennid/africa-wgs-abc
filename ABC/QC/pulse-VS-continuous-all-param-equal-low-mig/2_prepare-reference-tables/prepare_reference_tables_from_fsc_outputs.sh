#Gwenna Breton
#20220622
#Goal: once the fsc simulations have run, concatenate the different summary statistics to obtain a reference table where one row = sumstats from a simulation.
#This should handle cases when a simulation failed to provide sumstats (e.g. because the start parameters were not realistic).
#Caution! This script is to be taken as a guide and cannot be run per se.

folder=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/migration_models_Third_test/
cd ${folder}/output

###
### INITIAL RUNS
###

#Compute sumstats on ASD matrices.
L=199 #TODO adjust as needed
for scenario in scenario_1a_migration_asymmetric_0-0.0125 scenario_1b_migration_asymmetric_0-0.0125 scenario_1c_migration_asymmetric_0-0.0125 scenario_2a_migration_asymmetric_0-0.0125 scenario_2b_migration_asymmetric_0-0.0125 scenario_2c_migration_asymmetric_0-0.0125 scenario_3a_migration_asymmetric_0-0.0125 scenario_3b_migration_asymmetric_0-0.0125; do #TODO adjust as needed
for start in {1,201,401,601,801}; do #TODO adjust as needed
end=$(($start+$L))
for i in $(eval echo "{$start..$end}"); do
echo "tmp" > ${folder}/output/ASD/${scenario}_${i}_1_1.asd.dist.sumstats #That way if there is not ASD matrix, something gets written to the output and we have the right number of rows.
Rscript ${folder}/code/compute_statistics_based_on_ASD_forrackham.R ${scenario}_${i}_1_1 ${folder}/output/ASD/
paste ${folder}/output/ASD/${scenario}_${i}_1_1.asd.dist.sumstats >> ${folder}/output/ASDsumstats/${scenario}_ASDsumstats_${start}-${end}
rm ${folder}/output/ASD/${scenario}_${i}_1_1.asd.dist.sumstats
done
done
done

#Create reference tables: Merge to the param + sumstats files.
cd ${folder}/output/ref_tables

L=199
for scenario in scenario_1a_migration_asymmetric_0-0.0125 scenario_1b_migration_asymmetric_0-0.0125 scenario_1c_migration_asymmetric_0-0.0125 scenario_2a_migration_asymmetric_0-0.0125 scenario_2b_migration_asymmetric_0-0.0125 scenario_2c_migration_asymmetric_0-0.0125 scenario_3a_migration_asymmetric_0-0.0125 scenario_3b_migration_asymmetric_0-0.0125; do
for start in {1,201,401,601,801}; do
end=$(($start+$L))
paste ${folder}/output/sumstats/${scenario}_paramsets_sumstats_${start}-${end} ${folder}/output/ASDsumstats/${scenario}_ASDsumstats_${start}-${end} >> ${scenario}_paramsets_sumstats_ASDsumstats_1-1000
done
done

#Examine how many simulations have suceeded (i.e. simulations that do not have something like "File does not exist" or "tmp" in their vector of sumstats) (the "tmp" comes from the ASD matrices sumstats).
suffix=_paramsets_sumstats_ASDsumstats_
for range in 1-1000 ; do
for scenario in scenario_1a_migration_asymmetric_0-0.0125 scenario_1b_migration_asymmetric_0-0.0125 scenario_1c_migration_asymmetric_0-0.0125 scenario_2a_migration_asymmetric_0-0.0125 scenario_2b_migration_asymmetric_0-0.0125 scenario_2c_migration_asymmetric_0-0.0125 scenario_3a_migration_asymmetric_0-0.0125 scenario_3b_migration_asymmetric_0-0.0125; do
grep -v exist ${scenario}${suffix}${range} | grep -v tmp > ${scenario}${suffix}${range}_successful
done
done

#Get the indices for the simulations that were unsuccessful.
range=1-1000
suffix=_paramsets_sumstats_ASDsumstats_
for scenario in scenario_1a_migration_asymmetric_0-0.0125 scenario_1b_migration_asymmetric_0-0.0125 scenario_1c_migration_asymmetric_0-0.0125 scenario_2a_migration_asymmetric_0-0.0125 scenario_2b_migration_asymmetric_0-0.0125 scenario_2c_migration_asymmetric_0-0.0125 scenario_3a_migration_asymmetric_0-0.0125 scenario_3b_migration_asymmetric_0-0.0125; do
grep -e exist -e tmp ${scenario}${suffix}${range} > ${scenario}${suffix}${range}_unsuccessful
cut -f1 ${scenario}${suffix}${range}_unsuccessful | tr "\n" " " > ${scenario}${suffix}${range}_unsuccessful_iterations
done

#######################################
###
### Following this, new simulations are run for the vectors of parameters that have failed initially. In total, a given vector is run three times. If it fails the third time, a new vector is used instead.
### The outputs of the reruns are treated like above to identify the simulations that succeeded and the ones that failed. Then the sumstats for all successful iterations are concatenated.
### The files containing the simulations parameters and the sumstats have 421 columns, corresponding to:

#col1: iteration
#col2: seed
#col3-39: parameters of the model
#col40-421: sumstats

### These files are used to prepare the reference tables:

### 1. Reference table containing the sumstats:
##### Extract the sumstats:
cp /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/Second_trial/output/ref_tables_model_choice/Secondtrial_382sumstats_header.txt ABC_migration_thirdtest_8models_1000repeats_sumstats
for scenario in 1a 1b 1c 2a 2b 2c 3a 3b ; do
cut -f40-421 scenario_${scenario}_migration_asymmetric_0-0.0125${suffix}${range}_successful >> ABC_migration_thirdtest_8models_1000repeats_sumstats
done

##### Create a new vector where some statistics are proportions of bial.tot.nomiss. For simulated data, bial.tot.nomiss = bial.tot
R
data <- read.table(file="ABC_migration_thirdtest_8models_1000repeats_sumstats",header=TRUE)
prop_bial_hom = data[,3:12]/data[,1]
prop_private = data[,38:42]/data[,1]
write.table(cbind(prop_bial_hom,prop_private),file="ABC_migration_thirdtest_8models_1000repeats_sumstats_prop",col.names=c('prop.bial.NK','prop.bial.SK','prop.bial.WP','prop.bial.EP','prop.bial.NP','prop.hom.NK','prop.hom.SK','prop.hom.WP','prop.hom.EP','prop.hom.NP','prop.private.bial.NK','prop.private.bial.SK','prop.private.bial.WP','prop.private.bial.EP','prop.private.bial.NP'),quote=FALSE,row.names=FALSE)
q()
n

tr ' ' '\t' < ABC_migration_thirdtest_8models_1000repeats_sumstats_prop > tmp; mv tmp ABC_migration_thirdtest_8models_1000repeats_sumstats_prop

##### Put everything together:
root=ABC_migration_thirdtest_8models_1000repeats_sumstats
cut -f1-2 ${root} > col1-2
cut -f1-10 ${root}_prop > col3-12
cut -f13-37 ${root} > col13-37
cut -f11-15 ${root}_prop > col38-42
cut -f43-382 ${root} > col43-382
paste col1-2 col3-12 col13-37 col38-42 col43-382 > ABC_migration_thirdtest_8models_1000repeats_382sumstatswithprop
rm col*


### 2. File with model index.
i=1; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done #TODO adjust the value of i depending on the model index code, and the range for j depending on the number of repeats
i=2; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done
i=3; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done
i=4; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done
i=5; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done
i=6; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done
i=7; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done
i=8; for j in {1..1000}; do echo ${i} >> ABC_migration_thirdtest_8models_1000repeats_whichmodel; done


### 3. Files containing vectors of parameters:
#One file per model.
for scenario in 1a 1b 1c 2a 2b 2c 3a 3b ; do
cut -f3-39 scenario_${scenario}_migration_asymmetric_0-0.0125${suffix}${range}_3000repeats | cat ${folder}/models_and_parameters_sets/scenario_${scenario}_migration_asymmetric_0-0.0125.def_header - > scenario_${scenario}_migration_asymmetric_0-0.0125${suffix}${range}_3000repeats_params
done

### 4. File with correspondence between model number and model specifications
#This file is not used in any code, it is for personal reference. It should explain what the model code (e.g. "1") means.
#Example:
#1a	0-1	1
#1b	0-1	2
#1c	0-1	3
#2a	0-1	4
#2b	0-1	5
#2c	0-1	6
#3a	0-1	7
#3b	0-1	8
#Here model with code "1" corresponds to topology 1a and to migration rates taken in interval [0,1].

###
### Wrapup: Move everything to a folder, zip it.
###
zip -r reference_tables_continuous_0-0.0125_20210430_3000repeats.zip reference_tables_continuous_0-0.0125_20210430_3000repeats

