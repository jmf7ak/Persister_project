This folder contains the code for integrating the RNA-seq data with the metabolic models, performing gene essentiality and identifying genes for mutant experiment.

* persister_open_exchange.ipynb contains the code for processing the gene essentiality and flux sampling data to identify potential gene targets for persister cells.

* run_riptide_open.py and run_riptide_open_ATP.py contain the scripts to integrate the transcriptomics data into the iPau21 metabolic mode, generating flux samples and gene essentiality data.

*the submit_jobs.sh and submit_jobs_ATP.sh scripts  send all the riptide jobs to the UVA hpc for parallel processing.

*pyproject.toml contains the necessary requirements for persister_open_exchange.ipynb
*pyproject_riptide.toml contains the necessary requirements for the riptide codes on HPC
