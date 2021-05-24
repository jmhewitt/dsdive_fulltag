library(targets)

tar_make_future(
  names = c('validate_deep_surv_samples', 'dive_survival_eda_plot',
            'cee_dive_response_probs'), 
  workers = 400
)

tar_make(dive_survival_eda_plot, callr_function = NULL)
