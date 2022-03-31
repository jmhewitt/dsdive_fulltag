source('_targets.R')

sapply(as.list(tar_option_get('packages')), 
       function(x) library(x[[1]], character.only = TRUE))

tar_load(raw_sattags)
tar_load(data_pkg)
tar_load(deep_dive_depth)
tar_load(cape_hatteras_loc)
tar_load(covariate_tx_control)
shallow_threshold = 200
tar_load(template_bins)

# tar_load(fit_marginalized_model)

fit_marginalized_model = list(
  samples = file.path('output', 'mcmc', 'fit_marginalized_model'),
  package = file.path('output', 'mcmc', 'fit_marginalized_model', 'nim_pkg.rds')
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
burn = 1:(nrow(samples)*.1)

# get top-level names and groupings of variables being sampled
sampling_targets = colnames(samples)
sampling_target_groups = gsub(
  pattern = '\\[.*\\]',
  replacement = '',
  x = sampling_targets
)
sampling_groups = unique(sampling_target_groups)

#
# extract information about subject
#

# read in array parameter, which serves as target animal to test
subj_ind = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  
# wrap primary code in function
cee_preds = function() {
  
  # unwrap input
  raw_tag = raw_sattags[[subj_ind]]
  
  # get name for tag
  tagName = raw_tag$tag
  
  # determine subject id number for this tag
  subject_id = which(
    nim_pkg$consts$subject_id_labels == raw_tag$tag
  )
  
  # get info. for the subject's final segment that was analyzed
  final_seg = nim_pkg$consts$segments[
    max(which(nim_pkg$consts$segments[,'subject_id'] == subject_id)),
  ]
  
  # skip processing if the animal was not clearly exposed
  if(!is.finite(raw_tag$exposure_time)) {
    err_msg = paste('There is no exposure time recorded for',
                    tagName,
                    'so we cannot analyze CEE response')
    warning(err_msg)
    return(list(list(error = err_msg, tag = tagName)))
  }
  
  # skip processing if baseline period ends considerably before exposure
  if(as.numeric(raw_tag$exposure_time) - 
     data_pkg$data$times[final_seg['end_ind']] > nim_pkg$consts$tstep) {
    err_msg = paste('There is a large gap between the end of the baseline',
                    'data and exposure time for',
                    tagName,
                    'so we cannot analyze CEE response')
    warning(err_msg)
    return(list(list(error = err_msg, tag = tagName)))
  }
  
  # TODO: consider skipping processing, or modifying the processing, if we
  # can't get enough data to compute the lagged covariates
  
  # indices of the final segment
  seg_inds = final_seg['start_ind']:final_seg['end_ind']
  
  #
  # re-build model so that we can use it to generate transition matrices
  #
  
  # need to re-link compiled functions when running target in parallel
  source(file.path('R', 'util', 'statespace_tools', 'statespace_tools.R'))
  
  # uncompiled model
  mod = nimbleModel(
    code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data,
    inits = nim_pkg$inits, name = tar_name(), calculate = FALSE
  )
  
  # compile model
  cmod = compileNimble(mod)
  
  # verify model has a finite likelihood
  cmod$calculate()
  
  #
  # composition sample!
  #
  
  # identify posterior samples that will be used in the posterior analysis
  posterior_sample_inds = (1:nrow(samples))[-burn]
  # posterior_sample_inds = seq(
  #   from = min(posterior_sample_inds),
  #   to = max(posterior_sample_inds),
  #   length.out = 1e4
  # )
  
  # draw from posterior predictive distribution for time-till-deep
  baseline_deep_pred_samples = sapply(posterior_sample_inds, function(ind) {
    message(ind)
    
    # transfer posterior sample of model parameters to model object
    for(tgt_group in sampling_groups) {
      cmod[[tgt_group]] = samples[
        ind, sampling_targets[sampling_target_groups == tgt_group]
      ]
    }
    
    # update model components
    cmod$calculate()
    
    txmat_seq = stageTxMats(
      betas = cmod$beta_tx[final_seg['subject_id'],,,], 
      covariates = cmod$covariates[, seg_inds], 
      n_timepoints = length(seg_inds)
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
      subj = subject_id,
      covariate_tx_control = covariate_tx_control, 
      shallow_threshold = shallow_threshold
    )
    # extract times taken to reach a deep depth (sec)
    c(first_deep = min(which(
        template_bins$center[sim$depth_bins] > deep_dive_depth
      )) * nim_pkg$consts$tstep,
      second_deep = length(sim$depth_bins) * nim_pkg$consts$tstep)
  })
  
  # observed amount of time it took to see a deep depth post-exposure (sec)
  observed_time_to_deep = c(
    first_deep = min(which(
      raw_tag$depths[raw_tag$exposed == 1] > deep_dive_depth
    )) * nim_pkg$consts$tstep,
    second_deep = data.frame(
      deep = raw_tag$depths[raw_tag$exposed == 1] > deep_dive_depth,
      recovery =raw_tag$depths[raw_tag$exposed == 1] < shallow_threshold
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
    tag = tagName,
    deep_at_exposure = 
      data_pkg$data$covariates['depth', final_seg['end_ind']] >
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

f = file.path('output', 'cee', 'samples')
dir.create(path = f, showWarnings = FALSE, recursive = TRUE)

fname = paste('cee_predictions_raw_sattags_', subj_ind, '.rds', sep = '')

saveRDS(res, file = file.path(f, fname))
