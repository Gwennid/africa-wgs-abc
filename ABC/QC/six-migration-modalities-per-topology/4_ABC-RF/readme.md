# Model choice with ABC RF

[This R script](20210518_ABC_RF_scenario-1a_6migrationmodels.R) contains the different steps of the model choice procedure with ABC-RF, where the six models are of topology 1a. We performed two ABC-RF analyses: model choice among the six models; and model choice among the six models grouped in two groups (pulse versus continuous). The script is given as an example; the same analyses were performed for the other seven topologies.

The outputs of the script fall into two categories: QC (e.g. PCA of the summary statistics, goodness-of-fit, information relative to the ABC-RF procedure) and the model choice results: confusion matrices (see examples for [six models](confusionmatrix_6Models.txt) and [two groups](confusionmatrix_2Groups.txt)) and votes. The latter two outputs are plotted using [this script](20210520_plot_heatmaps_votes_ABCRFmodelchoice.R). The ABC-RF result plots are shown [here](../5_ABC-RF-result-plots). 
