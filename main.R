###########################################################
## Mixture and catalytic models fitted to simulated data ##
###########################################################


## load packages
library(ggplot2)
library(mixdist)
library(epitools)
library(serostat) # if using a version of R where this package cannot be downloaded,
                  # see spline.R file for the required function's source code
library(scales)
library(gridExtra)
library(grid)

## create output directory
dir.create("model_output")
dir.create("model_output/catalytic")
dir.create("model_output/mixture")


## Generation of simulated data ##
source("R/simulating_data_functions.R")
source("scripts/simulate_data.R")

## Calculate titre thresholds per simulated dataset ##
source("scripts/simulated_data_thresholds.R")

## to run the catalytic models
source("R/cat_model_functions.R") 
source("scripts/to_run_catalytic_models.R")

## to run the mixture model
source("R/mix_fitting_functions.R")
source("R/mix_calc_functions.R")
source("scripts/to_run_mixture_model.R") # N.B. likely will have to manually run through this script
                                   # (loop may throw error when trying to fit distributions to some 
                                   # of the data sets)

## model evaluation and comparison
source("scripts/explore_results_1.R") # mixture model
source("scripts/explore_results_2.R") # catalytic models
source("scripts/explore_results_3.R") # all models
