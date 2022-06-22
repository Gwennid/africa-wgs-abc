#Gwenna Breton
#Goal: Create reference tables to compare the models with continuous versus pulse migration. Here we want to do a QC and to compare the modalities with very low and fixed migration rates and with the rest of the parameters fixed. Based on a script from 20210521
#Example: 1000 iterations of topology 1a.

###
### Folder locations
###

continuous=/migration_models_6_QC_very-low/output/ref_tables/ #Replace
pulse=/pulse_models_6_QC_very-low/output/ref_tables/ #Replace


###
### Topology 1a
###

cd /ref_tables_continuousVSpulse/QC/constant_parameters #Replace (output folder)
cd 1a/

#Table1 (sumstats)
head -1001 ${pulse}ABC_pulse_constant_0.000001_1model_1000repeats_sumstats_382sumstatswithprop > ABC_scenario_1a_2models_1000repeats_sumstats_382sumstatswithprop
head -1001 ${continuous}ABC_continuous_constant_0.000001_model1a_1000repeats_382sumstatswithprop | tail -n 1000 >> ABC_scenario_1a_2models_1000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -1000 ${pulse}ABC_pulse_constant_0.000001_1model_1000repeats_whichmodel > ABC_scenario_1a_2models_1000repeats_whichmodel
head -1000 ${continuous}ABC_continuous_constant_0.000001_model1a_1000repeats_whichmodel >> ABC_scenario_1a_2models_1000repeats_whichmodel

#Table3 (key for model number)
#Prepared manually

zip -r reference_tables_scenario_1a_2models_constant_1000repeats.zip 1a
