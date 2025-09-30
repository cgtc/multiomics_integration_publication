# Exploring multiomics integration strategies

This repo is divided into two main directories:
  - 01_scripts: this contains scripts/helper functions for obtaining the different data types from AWS/GitHub.
  - notebooks: this contains the analysis notebooks for the different integration strategies.

The 00_data directory is not included in this repo due to the confidential nature of the data obtained from Plasticell's NK cell differentiation process.

## Scripts directory

There are 10 R scripts and 1 shell script in this directory, all for loading/processing LCMS, miRNA, Olink, and scRNA data in some way. The data_acquisition script is used within the integration analysis notebooks to load in the data using these.

They are included here to illustrate some of the processing steps that each modality underwent. Some processing steps are not documented in these scripts, but are described in the associated paper (e.g. scRNA/miRNA processing steps).

## Notebooks directory

The analysis notebooks (qmds and ipynbs) in this repository contain the application of five different multiomic integration strategies: 
  - Early (random forest model used on concatenated data) - early_integration.qmd
  - Intermediate (multiblock-sPLS from mixOmics) - intermediate_integration_mixomics.qmd
  - Intermediate (neural network using PyTorch) - intermediate_integration_neural_network.ipynb
  - Hierarchical (random forest model on pathway scores using PathIntegrate package) - hierachical_integration.ipynb
  - Late (LASSO regression on independent datasets) - late_integration.qmd

The rendered html's are also included in this directory.

The outputs from each of these notebooks is loaded by comparing_integration_outputs.qmd to compare features/model performance etc..

## Version control

The packages and package versions used for analysis in R and python can be found in the renv.lock and requirements.txt files respectively.
