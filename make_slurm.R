library(targets)

tar_make_future(
  names = c('cee_dive_response_probs'), 
  workers = 400
)

tar_make(dive_survival_eda_plot, callr_function = NULL)
