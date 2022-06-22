# Prepare reference tables

Once the simulations are completed, reference tables for the ABC analysis are prepared. This is done in several steps.

## Step 1: Prepare reference tables for each model

This is performed once the simulations have been run as described [here](../../../fsc-simulations/readme.md). The outputs of the fsc simulations are concatenated and certain columns are extracted -and in one case, modified- to prepare the reference tables necessary for ABC model choice and parameter estimation.

After this step, one has the following files:
- a file with 382 sumstats;
- a file saying which row comes from which model (only relevant for model choice);
- a file with the vectors of parameters (only used for parameter estimation);
- a file describing which model code corresponds to which topology, migration rate etc.

[This bash script](prepare_reference_tables_from_fsc_outputs.sh) shows how this can be done.

## Step 2: Prepare reference tables for ABC model choice

For example in this case, we are interested in comparing the pulse and the continuous migration for each topology. So we concatenate the files for pulse and continuous migration for a given topology. This is shown in [this script](prepare_reference_tables_modelchoice_pulseVScontinuous_QC_constant_parameters.sh).

The reference tables at this stage are found [here](../3_reference-tables/).

- [ ] Decide whether it is a good idea to include these reference tables, or whether I should not include them at all, or the tables used for model choice (if someome wants to reproduce model choice).

## Step 3: Modify files for use with the R abcrf package

An example is shown in [this script](20210521_Prepare_inputs_for_modelchoice_QC_constant_parameters.R). It gives two outputs, with respective suffixes ".SUMSTATALL" and ".MODINDEX".

- [ ] The tables are found on my work computer here: /Users/gwennabreton/Documents/Previous_work/PhD_work/P2_RHG_KS/ms/Supporting_information/QC_all_param_equal_very_low_mig/20210908_backup_modelchoice-modelchoice_pulseVScontinuous-QC-Constant_parameters/1a/
