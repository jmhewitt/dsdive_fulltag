# fit model on single node
targets::tar_make(
  names = fit_marginalized_model,
  callr_function = NULL
)

if(interactive()) {
  # identify trajectories to use for seeding posterior predictive simulations 
  # used to help interpret model parameters
  #   - must be run from an interactive session as it requires user validation
  targets::tar_make(
    names = parameter_interpretation_patterns,
    callr_function = NULL
  )
}

# prepare validation data
targets::tar_make(
  names = validation_config_nstep,
  callr_function = NULL
)

# # compute CEE probabilities in parallel
# targets::tar_make_future(
#   names = marginalized_model_cee_prediction,
#   workers = 400,
#   callr_function = NULL
# )