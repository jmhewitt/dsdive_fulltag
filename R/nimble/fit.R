fit = function(nim_pkg, nsamples, nthin, max_batch_iter = Inf, 
               max_batch_time = Inf, max_batch_mem = 1024^3,
               sample_dir) {
  # Build an MCMC sampler to approximate posterior inference.  Sampler will 
  # run in batches, periodically dumping samples to disk.  The time between 
  # dumps depends on the number of samples in each batch, and this is set to
  # balance memory and time costs.  Time can between batch dumps can be 
  # controlled indirectly by lowering max_batch_iter, and the size of each dump
  # can be controlled by lowering max_batch_mem.  The former is important for
  # guarding against loss of data caused by general system failures, while the 
  # latter is important for guarding against loss of data caused by 
  # out-of-memory crashes.
  #
  # Parameters:
  #  nim_pkg - data used to populate sampler
  #  nsamples - number of MCMC samples to return (after thinning)
  #  nthin - number of samples to thin
  #  max_batch_time - maximum time per batch (seconds), default is "no limit"
  #  max_batch_iter - maximum number of samples per batch
  #  max_batch_mem - maximum size of samples per batch (bytes), default is 1 GiB
  #  sample_dir - directory in which to save sample files
  # 
  # Return: 
  #  a list containing file names of sampler outputs
  
  # better initialization for latent stages
  ffbs_stages_c = compileNimble(ffbs_stages)
  for(i in 1:nim_pkg$consts$n_segments) {
    seg_start = nim_pkg$consts$segments[i, 'start_ind']
    seg_end = seg_start + nim_pkg$consts$segments[i, 'length'] - 1
    inds = seg_start:seg_end
    nim_pkg$data$stages[inds] = ffbs_stages_c(
      depths = nim_pkg$data$depths[inds], 
      n_timepoints = nim_pkg$consts$segments[i, 'length'], 
      transition_matrices = nim_pkg$data$transition_matrices, 
      n_bins = nim_pkg$consts$n_bins,
      n_stages = length(nim_pkg$consts$movement_types),  
      stage_map = nim_pkg$consts$movement_types, 
      alpha = nim_pkg$inits$alpha_mu, beta = nim_pkg$inits$beta_mu, 
      covariates = nim_pkg$data$covariates[,inds], 
      pi_discretization = nim_pkg$consts$pi_discretization, 
      n_pi = nim_pkg$consts$n_pi,
      n_lambda = nim_pkg$consts$n_lambda, 
      lambda_discretization = nim_pkg$consts$lambda_discretization, 
      betas_tx = nim_pkg$inits$betas_tx[,,1], 
      stage_supports = nim_pkg$data$stage_supports[,inds]
    )
  }
  
  # define and compile model
  model = nimbleModel(code = modelCode, constants = nim_pkg$consts, 
                      data = nim_pkg$data, inits = nim_pkg$inits, 
                      name = tar_name())
  model_c = compileNimble(model, projectName = tar_name())
  
  # initialize MCMC config
  conf = configureMCMC(model)
  
  # remove all samplers for forage period descent parameters
  stage_ind = which(
    colnames(nim_pkg$inits$alpha_mu) == 'deep_forage'
  )
  for(i in 1:nim_pkg$consts$n_covariates) {
    conf$removeSampler(
      paste('alpha_mu[', i, ', ', stage_ind, ']', sep ='')
    )
    conf$removeSampler(
      paste('alpha_var[', i, ', ', stage_ind, ']', sep ='')
    )
    for(j in 1:nim_pkg$consts$n_subjects) {
      conf$removeSampler(
        paste('alpha[', i, ', ', stage_ind, ', ', j, ']', sep ='')
      )
    }
  }
  
  # remove samplers for movement parameter effects not being estimated
  covariate_inds = which(rownames(nim_pkg$data$covariates) != 'intercept')
  for(i in covariate_inds) {
    for(j in 1:nim_pkg$consts$n_stages) {
      conf$removeSampler(paste('alpha_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('beta_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('alpha_var[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('beta_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        conf$removeSampler(paste('alpha[', i, ', ', j, ', ', k, ']', sep = ''))
        conf$removeSampler(paste('beta[', i, ', ', j, ', ', k, ']', sep = ''))
      }
    }
  }
  
  # remove samplers for stage transition coefficients not being estimated
  covariate_inds = which(rownames(nim_pkg$data$covariates) != 'intercept')
  for(i in covariate_inds) {
    for(j in 1:nim_pkg$consts$n_stage_txs) {
      conf$removeSampler(paste('betas_tx_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('betas_tx_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        conf$removeSampler(
          paste('betas_tx[', i, ', ', j, ', ', k, ']', sep = '')
        )
      }
    }
  }
  
  # latent stage vector samplers
  conf$removeSampler('stages')
  for(seg_ind in 1:nim_pkg$consts$n_segments) {
    start_ind = nim_pkg$consts$segments[seg_ind, 'start_ind']
    end_ind = start_ind + nim_pkg$consts$segments[seg_ind, 'length'] - 1
    conf$addSampler(
      target = paste('stages[', start_ind, ':', end_ind, ']', sep = ''),
      type = 'Stage'
    )
  }
  
  # set monitors, to partition output
  conf$resetMonitors()
  conf$addMonitors(c('betas_tx_mu', 'betas_tx_var', 'alpha_mu', 'beta_mu', 
                     'alpha_var', 'beta_var'))
  conf$addMonitors2('stages')
  
  # build and compile sampler
  mcmc = buildMCMC(conf)
  mcmc_c = compileNimble(mcmc, resetFunctions = TRUE, showCompilerOutput = TRUE,
                         projectName = tar_name())
  
  
  # size of sampler output for model parameters
  param_samples = as.matrix(mcmc_c$mvSamples)
  param_samples_labels = object_size(colnames(param_samples))
  param_samples_content = object_size(as.numeric(param_samples))
  
  # size of sampler output for latent stages
  stage_samples = as.matrix(mcmc_c$mvSamples2)
  stage_samples_labels = object_size(colnames(stage_samples))
  stage_samples_content = object_size(as.numeric(stage_samples))
  
  # maximum batch size to meet max_batch_mem constraints
  labels_overhead = param_samples_labels + stage_samples_labels
  content_per_sample = param_samples_content + stage_samples_content
  max_iter_by_mem = 
    # max samples that can be in batch wrt memory constraints
    floor(as.numeric(
      (max_batch_mem - labels_overhead) / content_per_sample
    )) * 
    # number of MCMC iterations required for each sample in output
    nthin
  
  # estimate maximum batch size to meet max_batch_time constraints
  if(is.finite(max_batch_time)) {
    # crude est. of time per iteration; can be wrong if samplers are adaptive
    tick = proc.time()[3]
    mcmc_c$run(niter = 1)
    time_per_iter = proc.time()[3] - tick
    max_iter_by_time = floor(max_batch_time / time_per_iter)
  } else {
    max_iter_by_time = Inf
  }
  
  # set iterations per batch to the minimum of the constraints
  batch_iter = min(max_batch_iter, max_iter_by_mem, max_iter_by_time)
  
  # create output directory
  out_dir = file.path(sample_dir, tar_name())
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # aggregate files generated
  sampler_files = list(
    labels = file.path(out_dir,
      paste(tar_name(), '_parameter_samples_column_labels.rds', sep = '')
    ),
    parameter_samples = c(),
    stage_samples = c()
  )
  
  # save labels
  saveRDS(colnames(param_samples), file = sampler_files$labels)
  
  # sample!
  total_it = nsamples * nthin
  remaining_it = total_it
  batch_num = 1
  while(remaining_it > 0) {
    # determine number of iterations to run
    batch_iter = min(remaining_it, batch_size)
    # run sampler
    mcmc_c$run(niter = batch_iter, reset = FALSE, resetMV = TRUE)
    # save parameter samples, reducing file size by removing labels
    param_samples = as.matrix(mcmc_c$mvSamples)
    attr(param_samples, 'dimnames') = NULL
    f_params = file.path(out_dir, 
      paste(tar_name(), '_parameter_samples_', batch_num, '.rds', sep ='')
    )
    saveRDS(param_samples, file = f_params)
    sampler_files$parameter_samples = c(sampler_files$parameter_samples, f_params)
    # save stage labels, reducing file by downsizing to integers
    stage_samples = as.matrix(mcmc_c$mvSamples2)
    stage_samples = matrix(as.integer(stage_samples), nrow = nrow(stage_samples))
    f_stages = file.path(out_dir,
      paste(tar_name(), '_stage_samples_', batch_num, '.rds', sep ='')
    )
    saveRDS(stage_samples, file = f_stages)
    sampler_files$stage_samples = c(sampler_files$stage_samples, f_stages)
    # update progress
    message(paste(
      Sys.time(), ' :: Sampling ', 
      round((total_it - remaining_it) / total_it * 100, 2), 
      '% complete', sep = ''
    ))
    # increment counters
    remaining_it = remaining_it - batch_iter
    batch_num = batch_num + 1
  }
  
  # return list of output files
  sampler_files
}