library(targets)

tar_make_future(
  names = c('validate_deep_surv_samples', 'dive_survival_eda_plot'), 
  workers = 400
)