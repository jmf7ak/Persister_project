This folder contains the code for integrating the RNA-seq data with the metabolic models, then perfroming gene/reaction essentiality and machine learning analyses.

the persister_transcriptomics_no_metab codes integrate the counts data with the metabolic models and perform gene/ reaction essentality. 
  -persister_transcriptomics_no_metab optimizes for biomass
  -persister_transcriptomics_no_metab optimizes for ATP
  
  Essentiality_heatmap codes take the gene/reaction essentiality and flux sampling data from the models to make figure 5 as well as set up the data to train the random forest classifier (classification learner app)
