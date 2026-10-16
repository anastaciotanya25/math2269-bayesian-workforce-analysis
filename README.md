# Bayesian Workforce Analysis

This repository contains the analysis for a Bayesian statistics project examining which workforce characteristics are associated with annual employment growth across Australian occupations.

The analysis uses Jobs and Skills Australia ANZSCO occupation data and fits a robust Bayesian multiple regression model in R and JAGS.

## Repository Structure

- `scripts/` – data preparation, Bayesian model specification, and model execution
- `data/` – processed datasets used for modelling
- `figures/` – descriptive, model, and posterior figures
- `runs/` – outputs from MCMC runs and prior sensitivity analysis

## Notes

The original Jobs and Skills Australia Excel workbook is not included in this repository.

The `DBDA2E-utilities.R` file used for Bayesian diagnostics and plotting is also not included and should be obtained separately from the course materials.