## Simulating genetic datasets for ABC

We simulated genetic datasets with fastsimcoal2 (Excoffier and Foll, [2011](https://doi.org/10.1093/bioinformatics/btr124), Excoffier et al., [2013](https://doi.org/10.1371/journal.pgen.1003905)). The genetic model was based on the model in Jay et al. ([2019](https://doi.org/10.1093/molbev/msz038)).

The models are described in the fastsimcoal files with extensions ".tpl" and ".est". They can be found in the corresponding subfolders of the `continuous` and `pulse` folders. Note that we did not use the ".est" files to generate vectors of parameters, as the models were too complex for it to be convenient. Thus vectors of parameters were generated with a custom Python script using the prior information in the ".est". These vectors are in files with extension ".def", found in the corresponding subfolder (one subfolder per migration intensity: `low`, `intermediate` and `high`). For the migration models with pulse and possibility of high migration intensity, we performed model choice and parameter estimation; the ".def" are found in subfolders of `pulse/def/high/`.

The outputs from fastsimcoal, in ".arp" format, were converted to plink TPED format and to a set of Python objects in order to calculate summary statistics. This is done with script `code/convert_fsc_output_to_tped_and_compute_sumstats_Feb2020.py`.

Summary statistics were calculated with plink (Purcell et al., [2007](https://doi.org/10.1086/519795)) v1.90b4.9, R v3.6.1, Python v2.7.15, bash, awk, scripts by Jay et al., [2019](https://doi.org/10.1093/molbev/msz038) (accessible [here](https://gitlab.inria.fr/ml_genetics/public/demoseq/-/tree/master)), and the [asd software](https://github.com/szpiech/asd). The scripts for calculating summary statistics are `code/convert_fsc_output_to_tped_and_compute_sumstats_Feb2020.py`, `code/summary_statistics.py` and `code/compute_statistics_based_on_ASD_forrackham.R`. The first two scripts are called from the script that performs the simulations, while the latter script is executed separately.

The full code needed to perform one or several simulations is present in `code/ABC_simul_modelchoice_April2021_pipeline_top.sh` and `code/ABC_simul_modelchoice_Feb2020_pipeline_bottom_20200514_USINGSCRATCH.sh`. The three placeholders `SCENARIO`, `START` and `TO` should be replaced (using sed for example). An example of values is: `scenario_1a_migration_asymmetric_0-0.001`, `1` and `500`.
In the working directory, the subfolder `code` should be copied. A directory called `models_and_parameters_sets` should be created, containing the following three files for each scenario of choice: `_raw.tpl`, `.def` and `.def_header`. Preferably the prefix of these three files should be identical. Moreover, the following directories should be created to contain the outputs: `output/ASD`, `output/ASDsumstats`, `output/sumstats`

#TODO if not done in the ms: explain the strategy with the reruns.

## TODO

- [ ] Include: custom Python script to generate vectors of parameters (do I have it?)
- [x] Include: script to convert fsc output (or should it go somewhere else?
- [x] Example of script to perform the full simulation?
- [ ] What else?
- [ ] When is output/ref_tables used?
