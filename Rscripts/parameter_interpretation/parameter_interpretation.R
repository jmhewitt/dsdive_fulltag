#
# set random seed for task relative to entire job
#

# set seed for job
RNGkind("L'Ecuyer-CMRG")
set.seed(2022)

# set seed for task
taskId = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
s <- .Random.seed
for (i in 1:taskId) {
  s <- parallel::nextRNGStream(s)
}
.GlobalEnv$.Random.seed <- s

# 
# workspace configuration from targets workflow
#

source('_targets.R')

sapply(as.list(tar_option_get('packages')), 
       function(x) library(x[[1]], character.only = TRUE))

tar_load(cape_hatteras_loc)
tar_load(covariate_tx_control)
shallow_threshold = 200
tar_load(template_bins)

# enrich covariate control list with lon/lat, for celestial calculations
covariate_tx_control$lon = cape_hatteras_loc['lon']
covariate_tx_control$lat = cape_hatteras_loc['lat']

#
# task details
#

# load dive prototypes
ptypes = readRDS(file.path('output', 'parameter_interpretation', 
                           'parameter_interpretation_patterns.rds'))

# define times from which posterior predictive sampling will begin; only the 
# day/night/moonlit properties are important
init_times = list(
  daytime = seq(
    from = strptime(x = '2020-06-16 0600', format = '%Y-%m-%d %H%M', 
                    tz = 'US/Eastern'),
    by = covariate_tx_control$obs_freq,
    length.out = ncol(ptypes)
  ),
  night_dark = seq(
    from = strptime(x = '2020-01-25 1800', format = '%Y-%m-%d %H%M', 
                    tz = 'US/Eastern'),
    by = covariate_tx_control$obs_freq,
    length.out = ncol(ptypes)
  ),
  night_moonlit = seq(
    from = strptime(x = '2020-02-08 1805', format = '%Y-%m-%d %H%M', 
                    tz = 'US/Eastern'),
    by = covariate_tx_control$obs_freq,
    length.out = ncol(ptypes)
  )
)

# configurations for posterior predictive sampling
config_inds = expand.grid(
  ptype = 1:nrow(ptypes), 
  time = 1:length(init_times),
  stage = 1:nrow(tar_read(movement_classes)$stage_defs)
)

# enumerate all of the posterior predictive sampling needs
manifest = expand.grid(
  fit_rep = 0:50,
  post_pred_task = 1:nrow(config_inds)
)


#
# load data and information needed for posterior predictive sampling
#

# tar_load(fit_marginalized_model)

sample_dir = file.path(
  'output', 'mcmc', paste('fit_marginalized_model_', manifest$fit_rep[taskId], 
                          sep = '')
)

fit_marginalized_model = list(
  samples = sample_dir,
  package = file.path(sample_dir, 'nim_pkg.rds')
)

#
# load, label, and merge posterior samples
#

mvSample_files = dir(
  path = fit_marginalized_model$samples, 
  pattern = 'mvSamples_[0-9]+', 
  full.names = TRUE
)

mvSample2_files = dir(
  path = fit_marginalized_model$samples, 
  pattern = 'mvSamples2_[0-9]+', 
  full.names = TRUE
)

samples = do.call(rbind, lapply(mvSample_files, readRDS))
samples2 = do.call(rbind, lapply(mvSample2_files, readRDS))

colnames(samples) = readRDS(dir(
  path = fit_marginalized_model$samples, 
  pattern = 'mvSamples_colnames',
  full.names = TRUE
))

colnames(samples2) = readRDS(dir(
  path = fit_marginalized_model$samples, 
  pattern = 'mvSamples2_colnames',
  full.names = TRUE
))

samples = cbind(samples, samples2)
rm(samples2)

nim_pkg = readRDS(fit_marginalized_model$package)

# set burn-in
burn = 1:(nrow(samples)*.5)

# get top-level names and groupings of variables being sampled
sampling_targets = colnames(samples)
sampling_target_groups = gsub(
  pattern = '\\[.*\\]',
  replacement = '',
  x = sampling_targets
)
sampling_groups = unique(sampling_target_groups)

#
# draw posterior predictive samples
#

# # test sampling start times to make sure celestial covariates are as expected
# covariate_tx(
#   covariates = rbind(
#     daytime = NA,
#     moonlit = NA,
#     depth = ptypes[1,],
#     time = init_times$night_moonlit
#   ), 
#   control = covariate_tx_control
# )


# wrap primary simulation code in function
#
# only summary data is returned, not full set of simulated trajectories
# 
# Parameters: 
#  ptype - prototype dive to use to initialize posterior predictive samples
#  times - times to associate with prototype dive observations
#  init_stage - character vector (i.e., 'slow_descent') specifying the movement
#    type from which sampling begins
pred_samples = function(ptype, times, init_stage) {
  
  # map prototype to template bins
  ptype_bins = sapply(ptype, function(d) {
    which.min(abs(d - template_bins$center))
  })
  
  #
  # re-build model so that we can use it to generate transition matrices
  #
  
  # uncompiled model
  mod = nimbleModel(
    code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data,
    inits = nim_pkg$inits, name = basename(tempfile('model')), calculate = FALSE
  )
  
  # compile model
  cmod = compileNimble(mod, showCompilerOutput = TRUE)
  
  # verify model has a finite likelihood
  cmod$calculate()
  
  #
  # composition sample!
  #
  
  # identify posterior samples that will be used in the posterior analysis
  posterior_sample_inds = (1:nrow(samples))[-burn]
  
  message('Sampling')
  
  # draw from posterior predictive distribution for time-till-deep
  baseline_deep_pred_samples = sapply(posterior_sample_inds, function(ind) {
    
    # transfer posterior sample of model parameters to model object
    for(tgt_group in sampling_groups) {
      cmod[[tgt_group]] = samples[
        ind, sampling_targets[sampling_target_groups == tgt_group]
      ]
    }
    
    #
    # simulation uses population-level parameters instead of random effects
    #
    
    # back transform parameters to natural scale
    cmod$depth_tx_mat = exp(cmod$depth_tx_mat)
    
    # simulate dive until a deep depth is reached
    sim = fwd_sim_dive(
      stage = which(rownames(nim_pkg$consts$stage_defs) == init_stage),
      depth_bins = ptype_bins,
      covariates = rbind(
        daytime = NA,
        moonlit = NA,
        depth = ptype,
        time = times
      ), 
      n_timepoints = 1e3,
      nim_pkg = nim_pkg,
      model = cmod,
      lon = cape_hatteras_loc['lon'],
      lat = cape_hatteras_loc['lat'],
      depth_threshold = covariate_tx_control$deep_depth,
      template_bins = template_bins,
      subj = -1,
      covariate_tx_control = covariate_tx_control, 
      shallow_threshold = shallow_threshold
    )
    # extract times taken to reach a deep depth (sec)
    c(first_deep = min(which(
        template_bins$center[sim$depth_bins] > covariate_tx_control$deep_depth
      )) * nim_pkg$consts$tstep,
      second_deep = length(sim$depth_bins) * nim_pkg$consts$tstep)
  })

  # package results
  list(list(
    ptype = ptype,
    times = times,
    init_stage = init_stage,
    baseline_deep_pred_samples = baseline_deep_pred_samples
  ))
}

# generate samples for task
cfg = config_inds[manifest$post_pred_task[taskId],]
res = pred_samples(
  ptype = ptypes[cfg$ptype,], 
  times = init_times[[cfg$time]], 
  init_stage = rownames(nim_pkg$consts$stage_defs)[cfg$stage]
)

# annotate configuration in output
res[[1]]$config = list(
  prototype = cfg$ptype,
  time = names(init_times)[cfg$time],
  init_stage = rownames(nim_pkg$consts$stage_defs)[cfg$stage]
)


#
# save output
#

f = file.path(
  'output', 'parameter_interpretation',
  paste('parameter_interpretation_', manifest$fit_rep[taskId], sep = ''), 
  'samples'
)

# make sure output directory is empty (i.e., clear previous model output)
unlink(x = fpath, recursive = TRUE)

# (re-)create output directory
dir.create(path = f, showWarnings = FALSE, recursive = TRUE)

fname = paste(
  'parameter_interpretation_samples_cfg', 
  manifest$post_pred_task[taskId], 
  '.rds', 
  sep = ''
)

saveRDS(res, file = file.path(f, fname))
