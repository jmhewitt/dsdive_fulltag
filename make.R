# # fit model on single node
# targets::tar_make(
#   names = fit_marginalized_model,
#   callr_function = NULL
# )

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