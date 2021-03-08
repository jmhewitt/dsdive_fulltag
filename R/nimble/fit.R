fit = function(nim_pkg, nsamples, nthin, max_batch_iter = Inf, 
               max_batch_time = Inf, max_batch_mem = 1024^3,
               sample_dir, structural_covariates) {
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
  #  structural_covariates - list of covariates for which effects should not 
  #    be estimated
  # 
  # Return: 
  #  a list containing file names of sampler outputs
  
  # define and compile model
  model = nimbleModel(code = modelCode, constants = nim_pkg$consts, 
                      data = nim_pkg$data, inits = nim_pkg$inits, 
                      name = tar_name())
  model_c = compileNimble(model, projectName = tar_name())
  
  # initialize MCMC config
  conf = configureMCMC(model)
  
  # remove all descent prob. samplers for stages w/fixed descent parameters
  stage_ind = which(
    colnames(nim_pkg$inits$alpha_mu) %in% c('deep_forage', 'free_surface')
  )
  for(i in 1:nim_pkg$consts$n_covariates) {
    conf$removeSampler(
      paste('alpha_mu[', i, ', ', stage_ind, ']', sep ='')
    )
    # conf$removeSampler(
    #   paste('alpha_var[', i, ', ', stage_ind, ']', sep ='')
    # )
    # for(j in 1:nim_pkg$consts$n_subjects) {
    #   conf$removeSampler(
    #     paste('alpha[', i, ', ', stage_ind, ', ', j, ']', sep ='')
    #   )
    # }
  }
  
  # do not allow most covariates to influence descent probs
  covariate_inds = which(
    !(rownames(nim_pkg$data$covariates) %in% c('intercept'))
  )
  for(i in covariate_inds) {
    for(j in 1:nim_pkg$consts$n_stages) {
      conf$removeSampler(paste('alpha_mu[', i, ', ', j, ']', sep = ''))
      # conf$removeSampler(paste('alpha_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        # conf$removeSampler(paste('alpha[', i, ', ', j, ', ', k, ']', sep = ''))
      }
    }
  }
  
  # do not estimate effects of static covariates on speed
  static_covariates = c(
    structural_covariates, 'deep_depth', 'shallow_depth', 
    'time_since_surface', 'depth'
  )
  covariate_inds = which(
    rownames(nim_pkg$data$covariates) %in% static_covariates
  )
  for(i in covariate_inds) {
    for(j in 1:nim_pkg$consts$n_stages) {
      conf$removeSampler(paste('beta_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('beta_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        conf$removeSampler(paste('beta[', i, ', ', j, ', ', k, ']', sep = ''))
      }
    }
  }
  
  # do not estimate effects of static covariates on stage tx params
  static_covariates = c(
    structural_covariates, 'deep_depth', 'shallow_depth', 
    'time_since_surface', 'depth'
  )
  covariate_inds = which(
    rownames(nim_pkg$data$covariates) %in% static_covariates
  )
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
  
  # identifiability constraint: depth doesn't influence surface-level stage tx.
  covariate_inds = which(
    rownames(nim_pkg$data$covariates) %in% c('depth')
  )
  stage_tx_inds = which(
    grepl(pattern = '(deep_ascent|shallow_ascent|free_surface)__', 
          x = colnames(nim_pkg$inits$betas_tx_mu))
  )
  for(i in covariate_inds) {
    for(j in stage_tx_inds) {
      # set variables to identity
      model_c[[paste('betas_tx_mu[', i, ', ', j, ']', sep = '')]] = 0
      model_c[[paste('betas_tx_var[', i, ', ', j, ']', sep = '')]] = 1
      # remove samplers
      conf$removeSampler(paste('betas_tx_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('betas_tx_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        # set variable to identity
        model_c[[paste('betas_tx[', i, ', ', j, ', ', k, ']', sep = '')]] = 0
        # remove sampler
        conf$removeSampler(
          paste('betas_tx[', i, ', ', j, ', ', k, ']', sep = '')
        )
      }
    }
  }
  
  # restricted stage tx model: don't estimate deprecated effects
  stage_tx_inds = which(
    colnames(nim_pkg$inits$betas_tx_mu) %in% 
      c('deep_ascent__deep_descent', 'shallow_ascent__deep_descent')
  )
  for(i in 1:nim_pkg$consts$n_covariates) {
    for(j in stage_tx_inds) {
      # set variables to identity
      model_c[[paste('betas_tx_mu[', i, ', ', j, ']', sep = '')]] = 0
      model_c[[paste('betas_tx_var[', i, ', ', j, ']', sep = '')]] = 1
      # remove samplers
      conf$removeSampler(paste('betas_tx_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('betas_tx_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        # set variables to identity
        model_c[[paste('betas_tx[', i, ', ', j, ', ', k, ']', sep = '')]] = 0
        # remove samplers
        conf$removeSampler(
          paste('betas_tx[', i, ', ', j, ', ', k, ']', sep = '')
        )
      }
    }
  }
  
  # # for identifiability, don't estimate deviations from stage tx offsets for 
  # # stages with only one viable transition.  this deconflicts population-level
  # # means while leaving random effects untouched.
  # stage_tx_inds = which(
  #   !(nim_pkg$consts$betas_tx_stage_from %in% nim_pkg$consts$ascent_like_stages)
  # )
  # for(i in 1:nim_pkg$consts$n_covariates) {
  #   for(j in stage_tx_inds) {
  #     conf$removeSampler(paste('betas_tx_mu[', i, ', ', j, ']', sep = ''))
  #   }
  # }
  
  # # use blocked samplers for stage transition effects between stages
  # for(stage in nim_pkg$consts$ascent_like_stages) {
  #   # stage transition effects associated with stage
  #   stage_tx_inds = which(nim_pkg$consts$betas_tx_stage_from == stage)
  #   # replace intercept random effect samplers with blocked samplers
  #   for(k in 1:nim_pkg$consts$n_subjects) {
  #     # identify nodes
  #     tgt = paste('betas_tx_eps[', nim_pkg$consts$intercept_covariate, ', ', 
  #                 stage_tx_inds, ', ', k, ']', sep = '')
  #     # remove default samplers
  #     conf$removeSamplers(tgt)
  #     # add block RW samplers
  #     conf$addSampler(target = tgt, type = 'RW_block')
  #   }
  # }
  
  # don't estimate population-level effects if so requested 
  if(nim_pkg$consts$population_effects == FALSE) {
    conf$removeSamplers(c('betas_tx_mu', 'betas_tx_var', 'beta_mu', 'beta_var'))
  }
  
  # latent stage vector samplers
  conf$removeSampler('stages')
  # for(seg_ind in 1:nim_pkg$consts$n_segments) {
  #   start_ind = nim_pkg$consts$segments[seg_ind, 'start_ind']
  #   end_ind = start_ind + nim_pkg$consts$segments[seg_ind, 'length'] - 1
  #   conf$addSampler(
  #     target = paste('stages[', start_ind, ':', end_ind, ']', sep = ''),
  #     type = 'Stage'
  #   )
  # }
  
  # reset monitors, to partition output
  conf$resetMonitors()
  # store stage samples separately from all other parameters
  conf$addMonitors2('stages')
  standard_params = setdiff(
    unique(gsub(
      pattern = "\\[.*\\]", replacement = '', 
      x = unlist(sapply(conf$getSamplers(), function(s) s$target))
    )),
    'stages'
  )
  conf$addMonitors(standard_params)
  
  # set thinning
  conf$setThin(nthin)
  conf$setThin2(nthin)
  
  # view final sampler configuration
  print(conf)
  
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
    active_samplers = file.path(out_dir,
      paste(tar_name(), '_active_samplers.rds', sep = '')
    ),
    pkg = file.path(out_dir, paste(tar_name(), '_nim_pkg.rds', sep = '')),
    parameter_samples = c(),
    stage_samples = c()
  )
  
  # save a copy of the input data
  saveRDS(nim_pkg, file = sampler_files$pkg)
  
  # save active sampler list
  saveRDS(sapply(conf$samplerConfs, function(s) s$target), 
          file = sampler_files$active_samplers)
  
  # save labels
  saveRDS(colnames(param_samples), file = sampler_files$labels)
  
  # sample!
  total_it = nsamples * nthin
  remaining_it = total_it
  batch_num = 1
  while(remaining_it > 0) {
    # determine number of iterations to run
    batch_chunk = min(remaining_it, batch_iter)
    # run sampler
    mcmc_c$run(niter = batch_chunk, reset = FALSE, resetMV = TRUE)
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
    # increment counters
    remaining_it = remaining_it - batch_chunk
    batch_num = batch_num + 1
    # update progress
    message(paste(
      Sys.time(), ' :: Sampling ', 
      round((total_it - remaining_it) / total_it * 100, 2), 
      '% complete', sep = ''
    ))
  }
  
  # return list of output files
  sampler_files
}