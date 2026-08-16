# Code for "Quantifying diets using Bayesian inference: a case study comparing juvenile fish diets among restored seascapes"

# Reproducing the analyses
 
This folder contains the scripts used to fit all models and generate all figures/tables reported in the manuscript and supplementary material. Scripts are numbered where execution order matters.

## Requirements
 
- R (>= 4.x) and the following packages: `dplyr`, `tidyr`, `rstan`, `ggplot2`, `purrr`, `readr`, `stringr`, `here`, `reshape2`, `PNWColors`, `ggpubr`, `ggmcmc`, `binom`, `scales`
- Stan model files in `stan/`: `logit_int.stan` (Frequency of Occurrence model), `GF_beta_reg.stan` (gut fullness model)
- Raw data at `data/rawData.csv`
- Helper functions sourced automatically from `R/` (plotting themes, `get_FO_diff()`, coefficient/PPC plotting functions)
All scripts use `here::here()` for paths, so run them from an R session opened at the repository root (or an RStudio project rooted there).

## 1. Frequency of Occurrence (FO) models
 
Run in order — **01 - 03 in the same R session**, as `03_run_logit_silverperch.R` combines results for the three species to build the pooled results table.
 
1. `01_run_logit_pinfish.R`
2. `02_run_logit_croaker.R`
3. `03_run_logit_silverperch.R`

## 2. Gut fullness model
 
Run `run_gut_fullness_proportional.R`. This fits the beta regression model for gut fullness.