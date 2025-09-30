# Exploring multiomics integration strategies

This repo is divided into two main directories:
  - 01_scripts: this contains scripts/helper functions for obtaining the different data types from AWS/GitHub.
  - notebooks: this contains the analysis notebooks for the different integration strategies.

## Notebooks directory

The analysis notebooks (qmds and ipynbs) in this repository contain the application of five different multiomic integration strategies: 
  - Early (random forest model used on concatenated data)
  - Intermediate (multiblock-sPLS from mixOmics)
  - Intermediate (neural network using PyTorch)
  - Hierarchical (random forest model on pathway scores using PathIntegrate package)
  - Late (LASSO regression on independent datasets)

## Version control

The packages and package versions used for analysis in R and python can be found in the renv.lock and requirements.txt files respectively.
