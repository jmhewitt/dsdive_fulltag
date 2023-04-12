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

tar_load(data_pkg)
tar_load(raw_sattags)
tar_load(cape_hatteras_loc)
tar_load(covariate_tx_control)
tar_load(deep_dive_depth)
shallow_threshold = 200
tar_load(template_bins)

# enrich covariate control list with lon/lat, for celestial calculations
covariate_tx_control$lon = cape_hatteras_loc['lon']
covariate_tx_control$lat = cape_hatteras_loc['lat']


#
# task details
#

# load information about tags
tar_load(tag_timelines)
# compile list of all exposures across all tags
cees = do.call(rbind, lapply(tag_timelines, function(x) {
  if(nrow(x$cee_segments) > 0) {
    cbind(tag = x$tag, x$cee_segments)
  }
}))
# extract the cee to be analyzed
cee = cees[taskId,]

#
# load data and information needed for posterior predictive sampling
#

# load burned-in posterior samples
samples = readRDS(
  file.path('output', 'mcmc', 'fixed_init_beta', 'cee_predictive_samples_stacked.rds')
)

nim_pkg = readRDS(
  file.path('output', 'mcmc', 'fixed_init_beta', 'fit_marginalized_model_1',
            'nim_pkg.rds')
)

# #
# # load data and information needed for posterior predictive sampling
# #
# 
# sample_dir = file.path('output', 'mcmc', 'fit_marginalized_model_0')
# 
# fit_marginalized_model = list(
#   samples = sample_dir,
#   package = file.path(sample_dir, 'nim_pkg.rds')
# )
# 
# #
# # load, label, and merge posterior samples
# #
# 
# mvSample_files = dir(
#   path = fit_marginalized_model$samples, 
#   pattern = 'mvSamples_[0-9]+', 
#   full.names = TRUE
# )
# 
# mvSample2_files = dir(
#   path = fit_marginalized_model$samples, 
#   pattern = 'mvSamples2_[0-9]+', 
#   full.names = TRUE
# )
# 
# samples = do.call(rbind, lapply(mvSample_files, readRDS))
# samples2 = do.call(rbind, lapply(mvSample2_files, readRDS))
# 
# colnames(samples) = readRDS(dir(
#   path = fit_marginalized_model$samples, 
#   pattern = 'mvSamples_colnames',
#   full.names = TRUE
# ))
# 
# colnames(samples2) = readRDS(dir(
#   path = fit_marginalized_model$samples, 
#   pattern = 'mvSamples2_colnames',
#   full.names = TRUE
# ))
# 
# samples = cbind(samples, samples2)
# rm(samples2)
# 
# nim_pkg = readRDS(fit_marginalized_model$package)
# 
# # set burn-in
# burn = 1:(nrow(samples)*.5)

# get top-level names and groupings of variables being sampled
sampling_targets = colnames(samples)
sampling_target_groups = gsub(
  pattern = '\\[.*\\]',
  replacement = '',
  x = sampling_targets
)
sampling_groups = unique(sampling_target_groups)

#
# analysis routine
#

# wrap primary code in analysis function
cee_preds = function() {

  # normalize tag name
  cee_tag = gsub(pattern = '\\_DUML', replacement = '', x = cee$tag)
  
  # extract raw tag data for task
  tag_map = sapply(raw_sattags, function(x) x$tag)
  raw_tag = raw_sattags[[which(tag_map == cee_tag)]]
  
  #
  # find cee segment within the model data package
  #

  # map subject name to id in data package
  subj_id = which(nim_pkg$consts$subject_id_labels == cee_tag)

  # validate time data appears consistent with data package
  if(length(nim_pkg$data$depth_bins) != length(data_pkg$data$times)) {
    stop('Cannot re-associate observation times with data package entries')
  }

  # augment segment information with start and end times
  segment_df = data.frame(
    nim_pkg$consts$segments,
    start_time = as.POSIXct(
      data_pkg$data$times[nim_pkg$consts$segments[, 'start_ind']],
      origin = '1970-01-01 00:00.00 UTC', 
      tz = 'UTC'
    ),
    end_time = as.POSIXct(
      data_pkg$data$times[nim_pkg$consts$segments[, 'end_ind']],
      origin = '1970-01-01 00:00.00 UTC', 
      tz = 'UTC'
    )
  )

  # identify the segment that preceeds the CEE
  pre_exposure_segment = data.frame(segment_df) %>% 
    mutate(seg_id = 1:n()) %>% 
    filter(
      # only look at this subject's segments
      subject_id == subj_id,
      # only consider a segment that ends immediately before the cee
      end_time <= cee$start,
      cee$start <= end_time + nim_pkg$consts$tstep
    ) %>% 
    select(seg_id) %>% 
    unlist()

  # early return b/c the CEE cannot be analyzed
  if(length(pre_exposure_segment) == 0) {
    return(
      list(list(
        tag = cee$tag,
        cee = cee$cee_id,
        error = 'Missing pre-exposure baseline data to initialize predictions.'
      ))
    )
  }

  # indices of the final pre-exposure segment
  seg_inds = seq(
    from = nim_pkg$consts$segments[pre_exposure_segment, 'start_ind'],
    to = nim_pkg$consts$segments[pre_exposure_segment, 'end_ind']
  )

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
  
  # need to re-link compiled functions when running target in parallel, but
  # should only be run after the model is compiled for portability
  source(file.path('R', 'util', 'statespace_tools', 'statespace_tools.R'))

  #
  # composition sample!
  #

  # identify posterior samples that will be used in the posterior analysis
  # posterior_sample_inds = (1:nrow(samples))[-burn]
  posterior_sample_inds = (1:nrow(samples))

  message('Sampling')

  # draw from posterior predictive distribution for time-till-deep
  baseline_deep_pred_samples = sapply(posterior_sample_inds, function(ind) {
    
    message(ind)

    # transfer posterior sample of model parameters to model object
    for(tgt_group in sampling_groups) {
      cmod[[tgt_group]] = samples[
        ind, sampling_targets[sampling_target_groups == tgt_group]
      ]
    }
    
    # return depth bin transition matrices to natural scale
    cmod$depth_tx_mat = exp(cmod$depth_tx_mat)

    txmat_seq = stageTxMats(
      betas = cmod$beta_tx[subj_id,,,], 
      covariates = cmod$covariates[, seg_inds], 
      n_timepoints = length(seg_inds), 
      log = FALSE
    )
    
    # posterior predictive distribution for latent state before exposure
    xf = finalPred2LayerPartialCompressed(
      obs_lik_dict = cmod$depth_tx_mat,
      obs = cmod$depth_bins[seg_inds] - 1,
      txmat_seq = txmat_seq,
      x0 = cmod$init_stages,
      num_obs_states = nim_pkg$consts$n_bins
    )

    # simulate dive until a deep depth is reached
    sim = fwd_sim_dive(
      stage = sample(
        x = 1:nim_pkg$consts$n_stages, size = 1, prob = xf
      ),
      depth_bins = nim_pkg$data$depth_bins[seg_inds],
      covariates = data_pkg$data$covariates[, seg_inds], 
      n_timepoints = 1e3,
      nim_pkg = nim_pkg,
      model = cmod,
      lon = cape_hatteras_loc['lon'],
      lat = cape_hatteras_loc['lat'],
      depth_threshold = deep_dive_depth,
      template_bins = template_bins,
      subj = subj_id,
      covariate_tx_control = covariate_tx_control, 
      shallow_threshold = shallow_threshold
    )
    # extract times taken to reach a deep depth (sec)
    c(first_deep = min(which(
        template_bins$center[sim$depth_bins] > deep_dive_depth
      )) * nim_pkg$consts$tstep,
      second_deep = length(sim$depth_bins) * nim_pkg$consts$tstep)
  })

  # analyze depth data associated with cee
  cee_behaviors = data.frame(
    depth = raw_tag$depths,
    time = raw_tag$times,
    gap_after = raw_tag$gap_after
  ) %>%  
    filter(
      # bracket data to cee window
      cee$start <= time,
      time <= cee$end
    ) %>% 
    mutate(
      before_gap = cumsum(gap_after) == 0
    ) %>%
    filter(
      # only analyze filtered observations consecutive in time
      before_gap == TRUE
    )
  
  # early return observations not available
  if(nrow(cee_behaviors) == 0) {
    return(
      list(list(
        tag = cee$tag,
        cee = cee$cee_id,
        error = 'Missing post-exposure data.'
      ))
    )
  }
  
  # observed amount of time it took to see a deep depth post-exposure (sec)
  observed_time_to_deep = c(
    first_deep = min(which(
      cee_behaviors$depth > deep_dive_depth
    )) * nim_pkg$consts$tstep,
    second_deep = data.frame(
      deep = cee_behaviors$depth > deep_dive_depth,
      recovery = cee_behaviors$depth < shallow_threshold
    ) %>% 
      mutate(ind = 1:n(),
             visited_deep = cumsum(deep) > 0,
             recovered_after_deep = cumsum(visited_deep * recovery) > 0) %>% 
      filter(deep & recovered_after_deep) %>% 
      select(ind) %>% 
      slice(1) %>%
      unlist() %>% 
      as.numeric() * nim_pkg$consts$tstep
  )
  
  # package results
  list(list(
    tag = cee$tag,
    cee = cee,
    deep_at_exposure = 
      data_pkg$data$covariates['depth', tail(seg_inds,1)] >
      deep_dive_depth,
    baseline_deep_pred_samples = baseline_deep_pred_samples,
    observed_time_to_deep = observed_time_to_deep,
    pval = c(
      first_deep_shorter = mean(
        baseline_deep_pred_samples['first_deep',] <= 
        observed_time_to_deep['first_deep']
      ),
      second_deep_shorter = mean(
        baseline_deep_pred_samples['second_deep',] <=
        observed_time_to_deep['second_deep']
      ),
      second_deep_longer = mean(
        observed_time_to_deep['second_deep'] <=
        baseline_deep_pred_samples['second_deep',]
      )
    )
  ))
}

# execute script
res = cee_preds()

#
# save output
#

f = file.path('output', 'cee', 'samples_stacked')

dir.create(path = f, showWarnings = FALSE, recursive = TRUE)

fname = paste('cee_predictions_', taskId, '.rds', sep = '')

saveRDS(res, file = file.path(f, fname))
