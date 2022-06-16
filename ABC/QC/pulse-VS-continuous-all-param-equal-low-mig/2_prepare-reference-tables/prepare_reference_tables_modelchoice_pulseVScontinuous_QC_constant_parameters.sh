#20210521
#Gwenna Breton
#Goal: Create reference tables to compare the models with continuous versus pulse migration. Here we want to do a QC and to compare the modalities with very low and fixed migration rates and with the rest of the parameters fixed.
#I am starting with 500 iterations of topology 1a.

###
### Folder locations
###

continuous=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/migration_models_6_QC_very-low/output/ref_tables/
pulse=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/pulse_models_6_QC_very-low/output/ref_tables/


###
### Topology 1a
###

cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/constant_parameters
cd 1a/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_1model_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_1a_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model1a_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_1a_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_1model_1000repeats_whichmodel > ABC_scenario_1a_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model1a_1000repeats_whichmodel >> ABC_scenario_1a_2models_1000repeats_whichmodel

#Table3 (key for model number)
#TODO

#TODO later: the parameter tables.


zip -r reference_tables_20210521_scenario_1a_2models_constant_1000repeats.zip 1000repeats


###
### Topology 1b
###

cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/constant_parameters
cd 1b/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_model1b_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_1b_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model1b_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_1b_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_model1b_1000repeats_whichmodel > ABC_scenario_1b_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model1b_1000repeats_whichmodel >> ABC_scenario_1b_2models_1000repeats_whichmodel

cd ../
zip -r reference_tables_20210524_scenario_1b_2models_constant_1000repeats.zip 1b


###
### Topology 1c
###

cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/constant_parameters/1c/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_model1c_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_1c_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model1c_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_1c_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_model1c_1000repeats_whichmodel > ABC_scenario_1c_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model1c_1000repeats_whichmodel >> ABC_scenario_1c_2models_1000repeats_whichmodel

cd ../
zip -r reference_tables_20210524_scenario_1c_2models_constant_1000repeats.zip 1c


###
### Topology 2a
###
cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/constant_parameters/2a/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_model2a_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_2a_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model2a_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_2a_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_model2a_1000repeats_whichmodel > ABC_scenario_2a_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model2a_1000repeats_whichmodel >> ABC_scenario_2a_2models_1000repeats_whichmodel

cd ../
zip -r reference_tables_20210526_scenario_2a_2models_constant_1000repeats.zip 2a


###
### Topology 2b
###
cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/constant_parameters/2b/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_model2b_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_2b_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model2b_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_2b_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_model2b_1000repeats_whichmodel > ABC_scenario_2b_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model2b_1000repeats_whichmodel >> ABC_scenario_2b_2models_1000repeats_whichmodel

cd ../
zip -r reference_tables_20210526_scenario_2b_2models_constant_1000repeats.zip 2b


###
### Topology 2c
###
cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/constant_parameters/2c/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_model2c_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_2c_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model2c_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_2c_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_model2c_1000repeats_whichmodel > ABC_scenario_2c_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model2c_1000repeats_whichmodel >> ABC_scenario_2c_2models_1000repeats_whichmodel

cd ../
zip -r reference_tables_20210526_scenario_2c_2models_constant_1000repeats.zip 2c


###
### Topology 3a
###
cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/constant_parameters/3a/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_model3a_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_3a_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model3a_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_3a_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_model3a_1000repeats_whichmodel > ABC_scenario_3a_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model3a_1000repeats_whichmodel >> ABC_scenario_3a_2models_1000repeats_whichmodel

cd ../
zip -r reference_tables_20210526_scenario_3a_2models_constant_1000repeats.zip 3a


#The same code was used for 3b, but I erased it.






























