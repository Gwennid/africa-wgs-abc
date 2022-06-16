#20210518
#Gwenna Breton
#Goal: Create reference tables to compare the models with continuous versus pulse migration. Here we want to do a QC and to examine the six modalities of migration (pulse/continuous * three priors) for each scenario in turn (i.e. eight scenario). I decided to extract the relevant rows from the final reference tables.

###
### Folder locations
###

continuous_low=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/migration_models_Fifth_test/output/ref_tables/reference_tables_continuous_0-0.001_20210517_5000repeats/
#${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop
#${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel
#${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber

pulse_low=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/pulse_models_Fifth_test/output/ref_tables/reference_tables_pulse_0-0.01_20210507_5000repeats/
#${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop
#${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel
#${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber

pulse_inter=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/Fourth_trial/output/ref_tables_model_choice/reference_tables_pulse_0-0.25_20210506_5000repeats/
#${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop
#${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel
#${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber

continuous_inter=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/migration_models_Third_test/output/ref_tables/reference_tables_continuous_0-0.0125_20210517_5000repeats/
#${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop
#${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel
#${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber

pulse_high=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/Second_trial/output/ref_tables_model_choice/reference_tables_pulse_0-1_20210506_5000repeats/
#${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop
#${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel
#${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber

continuous_high=/proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/migration_models_First_test/output/ref_tables/reference_tables_continuous_0-0.05_20210517_5000repeats/
#${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop
#${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel
#${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber


###
###
###

cd /proj/snic2020-2-10/uppstore2017183/b2012165_nobackup/private/Seq_project_cont/ABC/ref_tables_continuousVSpulse/QC/scenario_by_scenario

cd 1a/

#Table1 (sumstats)
head -5001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_1a_6models_5000repeats_sumstats_382sumstatswithprop
head -5001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1a_6models_5000repeats_sumstats_382sumstatswithprop
head -5001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1a_6models_5000repeats_sumstats_382sumstatswithprop
head -5001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1a_6models_5000repeats_sumstats_382sumstatswithprop
head -5001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1a_6models_5000repeats_sumstats_382sumstatswithprop
head -5001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1a_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -5000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel > ABC_scenario_1a_6models_5000repeats_whichmodel
head -5000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel >> ABC_scenario_1a_6models_5000repeats_whichmodel
head -5000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel >> ABC_scenario_1a_6models_5000repeats_whichmodel
head -5000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel >> ABC_scenario_1a_6models_5000repeats_whichmodel
head -5000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel >> ABC_scenario_1a_6models_5000repeats_whichmodel
head -5000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel >> ABC_scenario_1a_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber > ABC_scenario_1a_6models_5000repeats_correspondence_modelnumber
head -1 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber >> ABC_scenario_1a_6models_5000repeats_correspondence_modelnumber
head -1 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber >> ABC_scenario_1a_6models_5000repeats_correspondence_modelnumber
head -1 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber >> ABC_scenario_1a_6models_5000repeats_correspondence_modelnumber
head -1 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber >> ABC_scenario_1a_6models_5000repeats_correspondence_modelnumber
head -1 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber >> ABC_scenario_1a_6models_5000repeats_correspondence_modelnumber

#TODO later: the parameter tables.

cd ../1b/

#Table1 (sumstats)
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_1b_6models_5000repeats_sumstats_382sumstatswithprop
head -10001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_sumstats_382sumstatswithprop
head -10001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_sumstats_382sumstatswithprop
head -10001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_sumstats_382sumstatswithprop
head -10001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_sumstats_382sumstatswithprop
head -10001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_sumstats_382sumstatswithprop
head -10001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -10000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel | tail -n 5000 > ABC_scenario_1b_6models_5000repeats_whichmodel
head -10000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_whichmodel
head -10000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_whichmodel
head -10000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_whichmodel
head -10000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_whichmodel
head -10000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1b_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -2 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber | tail -n 1 > ABC_scenario_1b_6models_5000repeats_correspondence_modelnumber
head -2 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1b_6models_5000repeats_correspondence_modelnumber
head -2 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1b_6models_5000repeats_correspondence_modelnumber
head -2 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1b_6models_5000repeats_correspondence_modelnumber
head -2 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1b_6models_5000repeats_correspondence_modelnumber
head -2 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1b_6models_5000repeats_correspondence_modelnumber

cd ../1c/

#Table1 (sumstats)
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_1c_6models_5000repeats_sumstats_382sumstatswithprop
head -15001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_sumstats_382sumstatswithprop
head -15001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_sumstats_382sumstatswithprop
head -15001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_sumstats_382sumstatswithprop
head -15001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_sumstats_382sumstatswithprop
head -15001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_sumstats_382sumstatswithprop
head -15001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -15000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel | tail -n 5000 > ABC_scenario_1c_6models_5000repeats_whichmodel
head -15000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_whichmodel
head -15000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_whichmodel
head -15000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_whichmodel
head -15000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_whichmodel
head -15000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_1c_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -3 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber | tail -n 1 > ABC_scenario_1c_6models_5000repeats_correspondence_modelnumber
head -3 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1c_6models_5000repeats_correspondence_modelnumber
head -3 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1c_6models_5000repeats_correspondence_modelnumber
head -3 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1c_6models_5000repeats_correspondence_modelnumber
head -3 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1c_6models_5000repeats_correspondence_modelnumber
head -3 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_1c_6models_5000repeats_correspondence_modelnumber

cd ../2a/

#Table1 (sumstats)
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_2a_6models_5000repeats_sumstats_382sumstatswithprop
head -20001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_sumstats_382sumstatswithprop
head -20001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_sumstats_382sumstatswithprop
head -20001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_sumstats_382sumstatswithprop
head -20001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_sumstats_382sumstatswithprop
head -20001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_sumstats_382sumstatswithprop
head -20001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -20000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel | tail -n 5000 > ABC_scenario_2a_6models_5000repeats_whichmodel
head -20000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_whichmodel
head -20000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_whichmodel
head -20000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_whichmodel
head -20000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_whichmodel
head -20000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2a_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -4 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber | tail -n 1 > ABC_scenario_2a_6models_5000repeats_correspondence_modelnumber
head -4 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2a_6models_5000repeats_correspondence_modelnumber
head -4 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2a_6models_5000repeats_correspondence_modelnumber
head -4 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2a_6models_5000repeats_correspondence_modelnumber
head -4 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2a_6models_5000repeats_correspondence_modelnumber
head -4 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2a_6models_5000repeats_correspondence_modelnumber

cd ../2b/

#Table1 (sumstats)
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_2b_6models_5000repeats_sumstats_382sumstatswithprop
head -25001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_sumstats_382sumstatswithprop
head -25001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_sumstats_382sumstatswithprop
head -25001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_sumstats_382sumstatswithprop
head -25001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_sumstats_382sumstatswithprop
head -25001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_sumstats_382sumstatswithprop
head -25001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -25000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel | tail -n 5000 > ABC_scenario_2b_6models_5000repeats_whichmodel
head -25000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_whichmodel
head -25000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_whichmodel
head -25000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_whichmodel
head -25000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_whichmodel
head -25000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2b_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -5 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber | tail -n 1 > ABC_scenario_2b_6models_5000repeats_correspondence_modelnumber
head -5 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2b_6models_5000repeats_correspondence_modelnumber
head -5 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2b_6models_5000repeats_correspondence_modelnumber
head -5 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2b_6models_5000repeats_correspondence_modelnumber
head -5 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2b_6models_5000repeats_correspondence_modelnumber
head -5 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2b_6models_5000repeats_correspondence_modelnumber

cd ../2c/

#Table1 (sumstats)
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_2c_6models_5000repeats_sumstats_382sumstatswithprop
head -30001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_sumstats_382sumstatswithprop
head -30001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_sumstats_382sumstatswithprop
head -30001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_sumstats_382sumstatswithprop
head -30001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_sumstats_382sumstatswithprop
head -30001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_sumstats_382sumstatswithprop
head -30001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -30000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel | tail -n 5000 > ABC_scenario_2c_6models_5000repeats_whichmodel
head -30000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_whichmodel
head -30000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_whichmodel
head -30000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_whichmodel
head -30000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_whichmodel
head -30000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_2c_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -6 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber | tail -n 1 > ABC_scenario_2c_6models_5000repeats_correspondence_modelnumber
head -6 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2c_6models_5000repeats_correspondence_modelnumber
head -6 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2c_6models_5000repeats_correspondence_modelnumber
head -6 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2c_6models_5000repeats_correspondence_modelnumber
head -6 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2c_6models_5000repeats_correspondence_modelnumber
head -6 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_2c_6models_5000repeats_correspondence_modelnumber

cd ../3a/

#Table1 (sumstats)
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_3a_6models_5000repeats_sumstats_382sumstatswithprop
head -35001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_sumstats_382sumstatswithprop
head -35001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_sumstats_382sumstatswithprop
head -35001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_sumstats_382sumstatswithprop
head -35001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_sumstats_382sumstatswithprop
head -35001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_sumstats_382sumstatswithprop
head -35001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -35000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel | tail -n 5000 > ABC_scenario_3a_6models_5000repeats_whichmodel
head -35000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_whichmodel
head -35000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_whichmodel
head -35000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_whichmodel
head -35000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_whichmodel
head -35000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3a_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -7 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber | tail -n 1 > ABC_scenario_3a_6models_5000repeats_correspondence_modelnumber
head -7 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3a_6models_5000repeats_correspondence_modelnumber
head -7 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3a_6models_5000repeats_correspondence_modelnumber
head -7 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3a_6models_5000repeats_correspondence_modelnumber
head -7 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3a_6models_5000repeats_correspondence_modelnumber
head -7 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3a_6models_5000repeats_correspondence_modelnumber

cd ../3b/

#Table1 (sumstats)
head -1 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop > ABC_scenario_3b_6models_5000repeats_sumstats_382sumstatswithprop
head -40001 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_sumstats_382sumstatswithprop
head -40001 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_sumstats_382sumstatswithprop
head -40001 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_sumstats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_sumstats_382sumstatswithprop
head -40001 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_sumstats_382sumstatswithprop
head -40001 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_sumstats_382sumstatswithprop
head -40001 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_382sumstatswithprop | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_sumstats_382sumstatswithprop

#Table2 (model number)
head -40000 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_whichmodel | tail -n 5000 > ABC_scenario_3b_6models_5000repeats_whichmodel
head -40000 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_whichmodel
head -40000 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_whichmodel
head -40000 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_whichmodel
head -40000 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_whichmodel
head -40000 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_whichmodel | tail -n 5000 >> ABC_scenario_3b_6models_5000repeats_whichmodel

#Table3 (key for model number)
#Comment: it is not the best key since it does not say "pulse" or "continuous", but it has the model number and the priors differ so it is ok.
head -8 ${pulse_high}ABC_pulse_0-1_8models_5000repeats_correspondence_modelnumber | tail -n 1 > ABC_scenario_3b_6models_5000repeats_correspondence_modelnumber
head -8 ${pulse_inter}ABC_pulse_0-0.25_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3b_6models_5000repeats_correspondence_modelnumber
head -8 ${pulse_low}ABC_pulse_0-0.01_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3b_6models_5000repeats_correspondence_modelnumber
head -8 ${continuous_high}ABC_continuous_0-0.05_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3b_6models_5000repeats_correspondence_modelnumber
head -8 ${continuous_inter}ABC_continuous_0-0.0125_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3b_6models_5000repeats_correspondence_modelnumber
head -8 ${continuous_low}ABC_continuous_0-0.001_8models_5000repeats_correspondence_modelnumber | tail -n 1  >> ABC_scenario_3b_6models_5000repeats_correspondence_modelnumber


zip -r reference_tables_20210518_scenario_1a_6models_5000repeats.zip 1a
zip -r reference_tables_20210518_scenario_1b_6models_5000repeats.zip 1b
zip -r reference_tables_20210518_scenario_1c_6models_5000repeats.zip 1c
zip -r reference_tables_20210518_scenario_2a_6models_5000repeats.zip 2a
zip -r reference_tables_20210518_scenario_2b_6models_5000repeats.zip 2b
zip -r reference_tables_20210518_scenario_2c_6models_5000repeats.zip 2c
zip -r reference_tables_20210518_scenario_3a_6models_5000repeats.zip 3a
zip -r reference_tables_20210518_scenario_3b_6models_5000repeats.zip 3b

 



















