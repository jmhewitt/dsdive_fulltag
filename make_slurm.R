library(targets)

tar_make_future(
  names = c('validate_deep_surv_samples', 'validate_deep_surv_distn'), 
  workers = 400
)

