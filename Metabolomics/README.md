
`src/` folder contains the code for this project
`data/` folder contains the processed metabolomics data from Metabolon 
`analysis/` folder contains all intermediary files that are produced from `run_Metabolomics.R` that go into `run_Metabolomics_visualizations.R`

1. `run_Metabolomics.R` process the metabolomics data from Metabolon, and determines metabolites significantly changing across time and across samples
2. `run_Metabolomics_visualizations.R` takes the output from `run_Metabolomics.R` and creates figures
