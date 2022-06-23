# prepare data for model fitting: dependencies for fit_script.R
targets::tar_make(
  names = c('cape_hatteras_loc', 'covariate_tx_control', 'data_pkg', 
            'movement_classes')
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

# TODO: refactor parameter_interpretation script to allow for multiple starts

#
# TODO: transfer to a post-sampling make script
#

targets::tar_make(
  names = c('parameter_interpretation_plot', 'random_effect_plot'),
  callr_function = NULL
)

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