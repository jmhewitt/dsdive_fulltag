library(targets)

tar_make_future(
  names = c('nim_fit_val', 'nim_fit'), 
  workers = 400
)