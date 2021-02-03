fit = function(nim_pkg) {
  
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
  
  model = nimbleModel(code = modelCode, constants = nim_pkg$consts, 
                      data = nim_pkg$data, inits = nim_pkg$inits, 
                      name = 'fulltag')
  
  model_c = compileNimble(model, projectName = 'fulltag')
  
  conf = configureMCMC(model)
  
  # remove samplers for forage period descent parameters
  deep_forage_stage = which(
    names(nim_pkg$consts$movement_types) == 'deep_forage'
  )
  for(i in 1:nim_pkg$consts$n_covariates) {
    conf$removeSampler(
      paste('alpha_mu[', i, ', ', deep_forage_stage, ']', sep ='')
    )
    conf$removeSampler(
      paste('alpha_var[', i, ', ', deep_forage_stage, ']', sep ='')
    )
    for(j in 1:nim_pkg$consts$n_subjects) {
      conf$removeSampler(
        paste('alpha[', i, ', ', deep_forage_stage, ', ', j, ']', sep ='')
      )
    }
  }
  
  # remove samplers for movement parameter effects not being used
  for(i in 2:nim_pkg$consts$n_covariates) {
    for(j in 1:nim_pkg$consts$n_stages) {
      conf$removeSampler(paste('alpha_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('beta_var[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('beta_mu[', i, ', ', j, ']', sep = ''))
      conf$removeSampler(paste('alpha_var[', i, ', ', j, ']', sep = ''))
      for(k in 1:nim_pkg$consts$n_subjects) {
        conf$removeSampler(paste('alpha[', i, ', ', j, ', ', k, ']', sep = ''))
        conf$removeSampler(paste('beta[', i, ', ', j, ', ', k, ']', sep = ''))
      }
    }
  }
  
  # only the first row
  conf$removeSampler('betas_tx')

  # latent stage vector samplers
  conf$removeSampler('stages')
  for(seg_ind in 1:nim_pkg$consts$n_segments) {
    start_ind = nim_pkg$consts$segments[seg_ind, 'start_ind']
    end_ind = start_ind + nim_pkg$consts$segments[seg_ind, 'length'] - 1
    conf$addSampler(
      target = paste('stages[', start_ind, ':', end_ind, ']', sep = ''),
      type = 'Stage',
      control = list(
        betas_tx_node = paste('betas_tx[1:', nim_pkg$consts$n_covariates,
                              ', 1:', ncol(nim_pkg$inits$betas_tx), ', ',
                              nim_pkg$consts$segments[seg_ind, 'subject_id'],
                              ']', sep = ''),
        stage_supports = nim_pkg$data$stage_supports[, start_ind:end_ind]
      )
    )
  }
  
  conf$addMonitors('stages')
  
  mcmc = buildMCMC(conf)
  
  mcmc_c = compileNimble(mcmc, resetFunctions = TRUE, showCompilerOutput = TRUE,
                         projectName = 'fulltag')
  
  mcmc_c$run(niter = 10)
  samples = as.matrix(mcmc_c$mvSamples)
  
  
}