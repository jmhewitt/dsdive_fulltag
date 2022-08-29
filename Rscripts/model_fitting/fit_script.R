#
# set random seed for task relative to entire job
#

# set seed for job
RNGkind("L'Ecuyer-CMRG")
set.seed(2022)

# set seed for task
taskId = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
s <- .Random.seed
for (i in 0:taskId) {
  s <- parallel::nextRNGStream(s)
}
.GlobalEnv$.Random.seed <- s

# 
# workspace configuration from targets workflow
#

source('_targets.R')

sapply(as.list(tar_option_get('packages')), 
       function(x) library(x[[1]], character.only = TRUE))

#
# load data and information needed for posterior sampling
#

tar_load(cape_hatteras_loc)
tar_load(covariate_tx_control)

# enrich covariate control list with lon/lat, for celestial calculations
covariate_tx_control$lon = cape_hatteras_loc['lon']
covariate_tx_control$lat = cape_hatteras_loc['lat']

#
# do posterior sampling
#

source(file.path('Rscripts','model_fitting','fit_model.R'))

fit_model(
  data_pkg = tar_read(data_pkg), 
  covariate_tx_control = covariate_tx_control, 
  movement_classes = tar_read(movement_classes),
  random_inits = taskId > 0, # non-random starts are used for taskId == 0
  random_betas = FALSE, # start all covariate effects from 0 or I
  out_dir = file.path(
    'output', 'mcmc', 'fixed_beta_starts',
    paste('fit_marginalized_model_', taskId, sep = '')
  )
)
