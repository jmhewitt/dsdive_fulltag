library(targets)

# library(nimble)
# source('_targets.R')

# run workflow in current session, allowing for debugging breakpoints
targets::tar_make(dive_length_predictions, callr_function = NULL)

targets::tar_make(nim_fit, callr_function = NULL)

targets::tar_make(parameter_diagnostics, callr_function = NULL)

targets::tar_make_future(names = parameter_diagnostics, workers = 6)

targets::tar_make_future(names = posterior_imputation_plots, workers = 6)

targets::tar_make(names = posterior_imputation_plots, callr_function = NULL)

targets::tar_make(names = mle_speeds, callr_function = NULL)
