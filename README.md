# Mixture-catalytic-models

This repository contains code to recreate the simulation study in the paper 'Estimating dengue transmission intensity from serological data: a comparative analysis using mixture and catalytic models' (unpublished). 

## Repository overview
This repository contains scripts to:
- simulate age-stratified antibody titre (serology) datasets
- estimate antibody titre thresholds per dataset to classify individual titres as seropositive or seronegative
- fit a mixture model
- assess the mixture model ability to correctly specify parameters in serology datasets
- fit time-varying and time-constant force of infection (FOI) catalytic models
- compare and plot the uncertainty and bias in the FOI and seroprevalence estimates from each model

## Code
All analyses can be run from the main.R file, which calls scripts to run the analysis (contained in the scripts folder) and functions (contained in the R folder). The purpose of each R file is described below.

### Scripts
#### _simulate_data.R_
  
Simulate 200 serological datasets - each dataset is a linelist of individuals with age and simulated antibody titre level.
  
#### _simulated_data_thresholds.R_
  
Estimates optimal antibody titre thresholds per dataset to classify the titres as seropositive (indicating a previous infection) or seronegative (indictating no previous infection) which is necessary for the catalytic models.
  
#### _to_run_catalytic_models.R_

Main script to fit the catalytic models.

#### _plotting_catalytic_models_1.R_

Code to generate and save .csv files of the seroprevalence and FOI estimates, and to generate plots of the estimated age-specific seroprevalence.

#### _to_run_mixture_model.R_

Main script to fit the mixture model.

#### _plotting_mixture_model_1.R_

Code to generate plots of the mixture model fitted to the histogram of the simulated antibody titre distributions.

#### _plotting_mixture_model_2.R_

Code to generate and save .csv files of the seroprevalence and FOI estimates, and to generate plots of the estimated age-specific seroprevalence.

#### _explore_results_1.R_

Mixture model results. Assess the ability of the mixture model to correctly assign the distribution family (Weibull, normal, gamma) to the seronegative and seropositive compartments of the simulated datasets. Compare and plot the estimated mean titre scores, distribution standard deviations and the FOI to the 'true' simulated values. Further, calculate uncertainty and bias in the FOI and seroprevalence estimates.  

#### _explore_results_2.R_

Catalytic model results. Calculate and plot uncertainty and bias in the FOI and seroprevalence estimates.

#### _explore_results_3.R_

Mixture and Catalytic model output comparisons. Compare the estimated FOI and seroprevalence values, bias and uncertainty between the three models and generate comparison plots.


### R
#### _simulating_data_functions.R_

Code for the simulate() function to simulate the data and for the est_serostatus() function which is minimised to calculate classification thresholds in the _simulated_data_thresholds.R_ script

#### _cat_model_functions.R_

Functions to process the categorised data and fit the time-varying FOI and time-constant FOI catalytic models.

#### _mix_fitting_functions.R_
 
Functions for fitting the mixture model to the titre distributions (step 1 of mixture model).

#### _mix_calc_functions.R_

Functions to estimate seroprevalence and FOI from the estimated mixture model parameters (step 2 of mixture model).

#### _spline.R_

This file contains the code for the mpspline.fit function from the serostat package (reference: Kovac, T 2018)
