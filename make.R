# prepare data for model fitting: dependencies for fit_script.R
targets::tar_make(
  names = c('cape_hatteras_loc', 'covariate_tx_control', 'data_pkg', 
            'movement_classes')
)
