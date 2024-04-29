# Model choice with ABC RF

[This R script](20210521_ABC_RF_scenario-1a_2migrationmodels_constant_parameters.R) contains the different steps of the model choice procedure with ABC-RF, where the two models are of topology 1a, with pulse versus continuous migration. It is given as an example. The same analysis was performed for the other seven topologies. There are slight differences in which summary statistics are constant (and hence are removed from the analysis) for the different topologies (though the bulk is the same).

The outputs of the script fall into two categories: QC (e.g. PCA of the summary statistics, goodness-of-fit, information relative to the ABC-RF procedure) and the model choice results: confusion matrices (see an [example](confusionmatrix_2Groups.txt)) and votes. The latter two outputs are plotted using [this script](20210525_plot_heatmaps_votes.R). The ABC-RF result plots are shown [here](../5_ABC-RF-result-plots). 
