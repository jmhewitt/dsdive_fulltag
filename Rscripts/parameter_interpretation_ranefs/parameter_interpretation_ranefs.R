# parameter interpretation for random effects

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
# load data and information needed for posterior predictive sampling
#

samples = readRDS(
  file.path('output', 'mcmc', 'fixed_init_beta', 'cee_predictive_samples.rds')
)

nim_pkg = readRDS(
  file.path('output', 'mcmc', 'fixed_init_beta', 'fit_marginalized_model_1', 
            'nim_pkg.rds')
)

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
  time = which(
    names(init_times) == 'daytime'
  ),
  stage = which(
    rownames(tar_read(movement_classes)$stage_defs) == 'fast_descent'
  )
)

# enumerate all of the posterior predictive sampling needs
manifest = expand.grid(
  subject_id = 1:nim_pkg$consts$n_subjects,
  post_pred_task = 1:nrow(config_inds)
)

#
# draw posterior predictive samples
#

# get top-level names and groupings of variables being sampled
sampling_targets = colnames(samples)
sampling_target_groups = gsub(
  pattern = '\\[.*\\]',
  replacement = '',
  x = sampling_targets
)
sampling_groups = unique(sampling_target_groups)


# wrap primary simulation code in function
#
# only summary data is returned, not full set of simulated trajectories
# 
# Parameters: 
#  ptype - prototype dive to use to initialize posterior predictive samples
#  times - times to associate with prototype dive observations
#  init_stage - character vector (i.e., 'slow_descent') specifying the movement
#    type from which sampling begins
#  subject_id - subject number to indicate whose random effects to use
pred_samples = function(ptype, times, init_stage, subject_id) {
  
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
  # posterior_sample_inds = (1:nrow(samples))[-burn]
  posterior_sample_inds = (1:nrow(samples))
  
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
      subject_id = subject_id,
      covariate_tx_control = covariate_tx_control, 
      shallow_threshold = shallow_threshold
    )
    # extract times taken to reach a deep depth (sec)
    res = c(first_deep = min(which(
        template_bins$center[sim$depth_bins] > covariate_tx_control$deep_depth
      )) * nim_pkg$consts$tstep,
      second_deep = length(sim$depth_bins) * nim_pkg$consts$tstep)
    rm(sim)
    res
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
  init_stage = rownames(nim_pkg$consts$stage_defs)[cfg$stage], 
  subject_id = manifest$subject_id[taskId]
)

# annotate configuration in output
res[[1]]$config = list(
  prototype = cfg$ptype,
  time = names(init_times)[cfg$time],
  init_stage = rownames(nim_pkg$consts$stage_defs)[cfg$stage],
  subject_id = manifest$subject_id[taskId]
)


#
# save output
#

f = file.path(
  'output', 'parameter_interpretation', 'fixed_init_beta', 'random_effects', 
  'samples'
)

dir.create(path = f, showWarnings = FALSE, recursive = TRUE)

fname = paste(
  'parameter_interpretation_samples_task_', 
  taskId, 
  '.rds', 
  sep = ''
)

saveRDS(res, file = file.path(f, fname))

q(save = 'no', status = 0, runLast = FALSE)
