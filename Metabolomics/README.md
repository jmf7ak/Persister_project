Make sure the METABOLON_DATA_origscaled_R file is in your R working directory as well as all the codes in the SRC folder

1. run_Metabolomics.R to process the metabolomics data. Visualizes the raw data and determines metabolites significantly changing across time and across samples
2. run_MetabsForIntegration.R to determine which metabolites are unique to each condition and should be integrated into each model
3. run_Metabolomics_visualizations.R creates figures
4. integrated_metabs.R creates bar graph for metabolites integrated and not integrated into the models
