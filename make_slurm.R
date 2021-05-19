library(targets)

tar_make_future(
  names = c('nim_fit', 'validate_deep_surv_samples'), 
  workers = 400
)