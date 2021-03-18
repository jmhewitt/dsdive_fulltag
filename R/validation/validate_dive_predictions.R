validate_dive_predictions = function(post_output_dir, validation_dives, burn,
                                     n_timepoints) {
  # Parameters:
  #  n_timepoints - number of samples from dive to use in prediction
  
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
  prob_deep = lapply(validation_dives, function(d) {
    # skip processing if dive is too short
    if(length(d$depths) < n_timepoints) {
      return(NA)
    }
    # skip processing if we observe the dive to be deep
    if(any(d$covariates['deep_depth',inds] == 1)) {
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
    stage_preds = apply(param_samples[-burn,], 1, function(theta) {
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
        covariates = d$covariates[,inds], 
        pi_discretization = nim_pkg$consts$pi_discretization,
        n_pi = nim_pkg$consts$n_pi, 
        n_lambda =nim_pkg$consts$n_lambda,  
        lambda_discretization = nim_pkg$consts$lambda_discretization, 
        betas_tx = matrix(theta[betas_tx_nodes], 
                          nrow = nim_pkg$consts$n_covariates), 
        stage_supports = d$stage_supports[,inds],
        surface_bin = d$surface_bin[inds]
      )
      # return final prediction for stage
      tail(sampled_stages, 1)
    })
    # posterior predictive probability that dive is a deep dive
    mean(stage_preds %in% deep_stages)
  })
  
  # summarize findings
  res = cbind(
    n_timepoints = n_timepoints,
    prob_deep = as.numeric(prob_deep),
    true_deep = sapply(validation_dives, function(x) x$true_deep),
    tag = sapply(validation_dives, function(x) x$tag)
  )
  
  # dump backup
  saveRDS(res, file = file.path(
    post_output_dir, 
    paste(tar_name(), '.rds', sep = '')
  ))
  
  res
}
