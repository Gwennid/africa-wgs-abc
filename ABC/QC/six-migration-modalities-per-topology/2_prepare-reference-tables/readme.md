# Prepare reference tables

Once the simulations are completed, reference tables for the ABC analysis are prepared. This is done in several steps.

## Step 1: Prepare reference tables for each model

Cf the example for another QC [here](../pulse-VS-continuous-all-param-equal-low-mig/2_prepare-reference-tables/).

- [ ] Replace by link to the actual file once I have it!

## Step 2: Prepare reference tables for ABC model choice

For example in this case, we are interested in comparing the pulse and the continuous migration for each topology. So we concatenated the files for pulse and continuous migration for a given topology. We included 5000 repeats for each model. This is shown in [this script](prepare_reference_tables_modelchoice_pulseVScontinuous_QC_byscenario.sh).

The reference tables at this stage are found [here](../3_reference-tables/).

- [ ] Decide whether it is a good idea to include these reference tables, or whether I should not include them at all, or the tables used for model choice (if someome wants to reproduce model choice).
