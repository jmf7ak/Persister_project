For this code to work, you need to install the COBRA toolbox, Gurobi solver, the TIGER toolbox and OPTGPsampler for MATLAB. It is recommended to download the entire model folder ad set it as your current matlab directory.

persisterSampling_newmodel.m is the main code that will call all of the other functions. 
Run each section of the code indivdually. 

Outputs of persister_Sampling_newmodel
1. excel sheets with flux sampling data for both the ATP and biomass optimized models (for statistical analysis in R)
2. excel sheets with the sensitivity analysis on MADE parameters for the gene knockouts: the number of models in which each gene is essential

The sampling data from this code is then put into the code Sampling.R to run statiscal analysis 
