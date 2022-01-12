# fit model on single node
targets::tar_make(
  names = fit_marginalized_model,
  callr_function = NULL
)

# # compute CEE probabilities in parallel
# targets::tar_make_future(
#   names = marginalized_model_cee_prediction,
#   workers = 400,
#   callr_function = NULL
# )