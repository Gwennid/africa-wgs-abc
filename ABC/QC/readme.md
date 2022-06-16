# QC ABC

We performed checks at various stages of our ABC analyses. In this folder we gathered files relevant to testing different aspects of our ABC framework.
The files used to validate the main results of the article (for example test inherent to the ABC-RF model choice procedure) can be found in the folders corresponding to the analyses in question, one level up.

The following QC were performed:

## Does the set of genetic windows the summary statistics are computed on matter?

We compared three sets of windows.

- [ ] What kind of information do I want to give here? The coordinates of the windows?

## Do the pulse and continuous models give similar results when all non-migration parameters are identical and migration rates are very low?

[Code](pulse-VS-continuous-all-param-equal-low-mig)

- [ ] This QC was used as an "example" thus other folders might refer to the files in there.
- [ ] TODO write short background (though there is something in the MMs already)
- [ ] Should I include the ABC-RF-QC plots? (PCA, goodness-of-fit, error of the forest...) (In the mean time they can be found here on my work computer: /Users/gwennabreton/Documents/Previous_work/PhD_work/P2_RHG_KS/ms/Supporting_information/QC_all_param_equal_very_low_mig/20210908_backup_modelchoice-modelchoice_pulseVScontinuous-QC-Constant_parameters/)

## Do the results of model choice between pulse and continuous migration for each topology match the results across topologies?

[Code](six-migration-modalities-per-topology)

- [ ] TODO write short background (though there is something in the MMs already)

## Is the proportion of migrants comparable in the pulse and continuous migration models?
