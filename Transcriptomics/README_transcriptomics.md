The code in this folder is used to analyze the counts from the RNA-seq data, which can be cound in the GEO database here:

1. DESeq2.R reads the counts files and determines the differential expression for each gene
2. run_RE_prepintegration.R takes the output files from DESeq2, filters the data for metabolic ganes and outputs a file for transcriptomic integration into the model via MADE
3. run_DE_visualizations.R creates the graphics for the paper from the differential expression data
