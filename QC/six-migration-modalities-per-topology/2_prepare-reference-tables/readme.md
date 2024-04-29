# Prepare reference tables

Once the simulations are completed, reference tables for the ABC analysis are prepared. This is done in several steps.

## Step 1: Prepare reference tables for each model

Cf the example for another QC [here](../pulse-VS-continuous-all-param-equal-low-mig/2_prepare-reference-tables/).

- [ ] Replace by link to the actual file once I have it!

## Step 2: Prepare reference tables for ABC model choice

For example in this case, we are interested in comparing the pulse and the continuous migration for each topology. So we concatenated the files for pulse and continuous migration for a given topology. We included 5000 repeats for each model. This is shown in [this script](prepare_reference_tables_modelchoice_pulseVScontinuous_QC_byscenario.sh).

The reference tables at this stage are too large to be uploaded. They can be found on my work computer here: /Users/gwennabreton/Documents/Previous_work/PhD_work/P2_RHG_KS/ms/Supporting_information/QC_by_topo/20210908_backup_modelchoice-modelchoice_pulseVScontinuous-QC-Scenario_by_scenario

- [ ] Decide whether I should try to include the reference tables or whether I skip it (see above).

## Step 3: Modify files for use with the R abcrf package

Cf the example for another QC [here](../pulse-VS-continuous-all-param-equal-low-mig/2_prepare-reference-tables/20210521_Prepare_inputs_for_modelchoice_QC_constant_parameters.R).

Parameter `nbSim` was set to `5000` and parameter `nbModels` was set to `6`.
