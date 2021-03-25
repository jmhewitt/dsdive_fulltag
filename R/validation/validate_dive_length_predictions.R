validate_dive_length_predictions = function(
  post_output_dir, validation_dives, burn, n_timepoints, timestep, lon, lat,
  template_bins
) {
  # Parameters:
  #  n_timepoints - number of samples from dive to use in prediction
  
  if(!exists('ffbs_stages_c')) {
    ffbs_stages_c = nimble::compileNimble(ffbs_stages)
  }
  
  # define paths to MCMC output
  paths = list(
    # location of mcmc output files
    path_out = post_output_dir,
    # identifying characteristics of key output files
    active_samplers = '*active_samplers.rds',
    param_sample_pattern = '*parameter_samples_[0-9]+',
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
  
  # ids for deep dive stages
  deep_stages = which(grepl('deep', names(nim_pkg$consts$movement_types)))
   
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
  
  # number of points in dive to use for prediction
  inds = 1:n_timepoints
  
  # test each dive
  remaining_time_samples = lapply(validation_dives, function(d) {
    # skip processing if dive is too short
    if(length(d$depths) < n_timepoints) {
      return(NA)
    }
    # to reduce computational requirements: skip processing if dive was not deep
    if(d$true_deep == FALSE) {
      return(NA)
    }
    # get modeled tag id 
    tag_id = which(d$tag == nim_pkg$consts$subject_id_labels)
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
    remaining_duration_preds = apply(param_samples[-burn,], 1, function(theta) {
      # draw stages from posterior predictive distribution
      sampled_stages = ffbs_stages_c(
        depths = d$depths[inds], 
        n_timepoints = length(inds),
        transition_matrices = nim_pkg$data$transition_matrices,
        n_bins = nim_pkg$consts$n_bins, 
        n_stages = nim_pkg$consts$n_stages, 
        stage_map = nim_pkg$consts$movement_types, 
        alpha = matrix(theta[alpha_nodes], nrow = nim_pkg$consts$n_covariates), 
        beta = matrix(theta[beta_nodes], nrow = nim_pkg$consts$n_covariates),
        covariates = d$covariates[, inds, drop = FALSE], 
        pi_discretization = nim_pkg$consts$pi_discretization,
        n_pi = nim_pkg$consts$n_pi, 
        n_lambda =nim_pkg$consts$n_lambda,  
        lambda_discretization = nim_pkg$consts$lambda_discretization, 
        betas_tx = matrix(theta[betas_tx_nodes], 
                          nrow = nim_pkg$consts$n_covariates), 
        stage_supports = d$stage_supports[, inds, drop = FALSE],
        surface_bin = d$surface_bin[inds]
      )
      # start forward simulation from final prediction for stage
      last_ind = tail(inds, 1)
      pred = fwd_sim_dive(
        stages = tail(sampled_stages, 1), depths = d$depths[last_ind], 
        covariates = d$covariates[, last_ind, drop = FALSE], 
        n_timepoints = Inf, nim_pkg = nim_pkg, 
        alpha = matrix(theta[alpha_nodes], nrow = nim_pkg$consts$n_covariates), 
        beta = matrix(theta[beta_nodes], nrow = nim_pkg$consts$n_covariates), 
        betas_tx = matrix(theta[betas_tx_nodes], 
                          nrow = nim_pkg$consts$n_covariates), 
        template_bins = template_bins, 
        times = as.POSIXct(x = nim_pkg$data$times[last_ind], 
                           origin = '1970-01-01 00:00.00 UTC', tz = 'UTC'), 
        timestep = timestep, lon = lon, lat = lat
      )
      # return posterior predictive sample of remaining time in dive (sec)
      diff(as.numeric(range(pred$times)))
    })
    # posterior predictive draws of remaining time in dive
    remaining_duration_preds
  })
  
  true_times = sapply(validation_dives, function(x) {
    # number of points in dive to use for prediction
    inds = 1:n_timepoints
    # true remaining time in dive (sec)
    diff(as.numeric(range(x$times[-inds])))
  })
  
  # package posterior samples for later assessment
  res = list(list(
    n_timepoints = n_timepoints,
    remaining_time_samples = remaining_time_samples,
    tag = sapply(validation_dives, function(x) x$tag),
    true_remaining_time = true_times,
    dive_start_time = sapply(validation_dives, function(x) x$dive_start_time),
    true_deep = sapply(validation_dives, function(x) x$true_deep)
  ))
  
  # dump backup
  saveRDS(res, file = file.path(
    post_output_dir, 
    paste(tar_name(), '.rds', sep = '')
  ))
  
  res
}
