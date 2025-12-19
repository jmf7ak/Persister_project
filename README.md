# Persister_project

This repository contains all code and data required to reproduce the results and figures from:

"Integrating multi-omics data with network modeling to characterize Pseudomonas aeruginosa persister cell metabolism." by Blazier, Ficarrotta, Metris, Kolling and Papin.

---

## Repository Structure

The project is organized into the following folders:

- `Metabolomics/` – Scripts and data for metabolomics analyses.
- `Modeling/` – Code for building and analyzing contextualized metabolic models.
- `Mutant_Validation/` – Scripts for analyzing single-gene mutant experimental data.
- `Transcriptomics/` – Data processing and analyses for RNA-sequencing, including differential expression and term enrichment.

Each folder contains a **README** with detailed instructions for running the code specific to that module.


The R environment is managed with `renv` to ensure reproducibility. The `renv.lock` file contains all package versions used in the analyses.
There are two py.toml files in `Modeling/`, `pyproject_riptide.toml` for creating a metabolic modeling environment and `pyproject.toml` for data analysis.
