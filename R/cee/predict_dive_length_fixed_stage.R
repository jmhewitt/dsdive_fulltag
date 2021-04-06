predict_dive_length_fixed_stage = function(
  post_output_dir, tag, burn, timestep, lon, lat, template_bins, fixed_stage
) {
  
  # define paths to MCMC output
  paths = list(
    # location of mcmc output files
    path_out = post_output_dir,
    # identifying characteristics of key output files
    active_samplers = '*active_samplers.rds',
    param_sample_pattern = '*parameter_samples_[0-9]+',
    stage_sample_pattern = 'stage_samples_[0-9]+',
    param_sample_column_labels = '*parameter_samples_column_labels.rds'
  )
  
  # load input used to generate MCMC output
  nim_pkg = readRDS(dir(path = paths$path_out, pattern = 'pkg', 
                        full.names = TRUE))
  
  # find posterior parameter samples 
  param_sample_files = dir(path = paths$path_out, 
                           pattern = paths$param_sample_pattern,
                           full.names = TRUE)
  
  # load posterior parameter samples 
  param_samples = do.call(rbind, lapply(param_sample_files, function(f) {
    readRDS(f)
  }))
  
  # label posterior parameter samples 
  colnames(param_samples) = readRDS(
    dir(path = paths$path_out, pattern = paths$param_sample_column_labels, 
        full.names = TRUE)
  )
  
  # clean up
  rm(param_sample_files)
   
  # common identifiers for model parameters
  gr = expand.grid(1:nim_pkg$consts$n_covariates, 1:nim_pkg$consts$n_stages)
  gr2 = expand.grid(1:nim_pkg$consts$n_covariates, 
                    1:nim_pkg$consts$n_stage_txs)
  alpha_indices = lapply(
    1:nrow(gr), 
    function(i) { 
      paste(gr[i,], collapse = ', ')
    }
  )
  alpha_nodes = paste('alpha_mu[', alpha_indices, ']', sep = '')
  
  # get modeled tag id 
  tag_id = which(tag == nim_pkg$consts$subject_id_labels)
  
  # get index for last modeled observation
  last_seg_ind = max(which(nim_pkg$consts$segments[, 'subject_id'] == tag_id))
  last_obs_ind = nim_pkg$consts$segments[last_seg_ind, 'end_ind']
  
  # build local identifiers for model parameters
  beta_indices = lapply(
    1:nrow(gr), 
    function(i) { 
      paste(c(gr[i,], tag_id), collapse = ', ')
    }
  )
  beta_nodes = paste('beta[', beta_indices, ']', sep = '')
  betas_tx_indices = lapply(
    1:nrow(gr2), 
    function(i) { 
      paste(c(gr2[i,], tag_id), collapse = ', ')
    }
  )
  betas_tx_nodes = paste('betas_tx[', betas_tx_indices, ']', sep = '')
  
  # sample over posterior distribution
  remaining_duration_preds = sapply((1:nrow(param_samples))[-burn], function(ind) {
    theta  = param_samples[ind,]
    # start forward simulation from final observation, using fixed stage
    pred = fwd_sim_dive(
      stages = fixed_stage, depths = nim_pkg$data$depths[last_obs_ind], 
      covariates = nim_pkg$data$covariates[, last_obs_ind, drop = FALSE], 
      n_timepoints = Inf, nim_pkg = nim_pkg, 
      alpha = matrix(theta[alpha_nodes], nrow = nim_pkg$consts$n_covariates), 
      beta = matrix(theta[beta_nodes], nrow = nim_pkg$consts$n_covariates), 
      betas_tx = matrix(theta[betas_tx_nodes], 
                        nrow = nim_pkg$consts$n_covariates), 
      template_bins = template_bins, 
      times = as.POSIXct(x = nim_pkg$data$times[last_obs_ind], 
                         origin = '1970-01-01 00:00.00 UTC', tz = 'UTC'), 
      timestep = timestep, lon = lon, lat = lat
    )
    # return posterior predictive sample of remaining time in dive (sec)
    diff(as.numeric(range(pred$times)))
  })
  
  # true remaining time
  
  # package posterior samples for later assessment
  res = list(list(
    remaining_time_samples = remaining_duration_preds,
    tag = tag, 
    pre_cee_time = as.POSIXct(x = nim_pkg$data$times[last_obs_ind], 
                              origin = '1970-01-01 00:00.00 UTC', tz = 'UTC')
  ))
  
  # dump backup
  saveRDS(res, file = file.path(
    post_output_dir, 
    paste(tar_name(), '.rds', sep = '')
  ))
  
  res
}
