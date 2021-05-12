nimble_targets = list(
  
  # define discrete movement classes/states, in which each class is associated
  # with one of the modeled lambda values, and one of the discretized pi values
  tar_target(
    name = movement_classes, 
    command = {
        res = list(
          lambda_class = c('slow', 'medium', 'fast'),
          pi_ind = c(
            'descent' = which(names(parameter_discretization$pi) == 'descent'),
            'free' = which(names(parameter_discretization$pi) == 'free'),
            'ascent' = which(names(parameter_discretization$pi) == 'ascent')
          )
      )
        
      res$stage_defs = rbind(
        c(pi_ind = as.numeric(res$pi_ind['descent']), 
          lambda_class = which(res$lambda_class == 'slow')),
        c(pi_ind = as.numeric(res$pi_ind['descent']), 
          lambda_class = which(res$lambda_class == 'medium')),
        c(pi_ind = as.numeric(res$pi_ind['descent']), 
          lambda_class = which(res$lambda_class == 'fast')),
        c(pi_ind = as.numeric(res$pi_ind['free']), 
          lambda_class = which(res$lambda_class == 'slow')),
        c(pi_ind = as.numeric(res$pi_ind['ascent']), 
          lambda_class = which(res$lambda_class == 'slow')),
        c(pi_ind = as.numeric(res$pi_ind['ascent']), 
          lambda_class = which(res$lambda_class == 'medium')),
        c(pi_ind = as.numeric(res$pi_ind['ascent']), 
          lambda_class = which(res$lambda_class == 'fast'))
      )
      
      rownames(res$stage_defs) = c(
        'slow_descent', 'medium_descent', 'fast_descent', 'slow_free', 
        'slow_ascent', 'medium_ascent', 'fast_ascent'
      )
        
      res
    }
  ),
  
  tar_target(
    name = nim_pkg, 
    command = flatten_tags(
      template_bins = template_bins, 
      lambda_discretization = parameter_discretization$lambda, 
      stage_defs = movement_classes$stage_defs, 
      init_movement_coefficients = list(lambda = rep(1,3)), 
      transition_matrices = transition_matrices, 
      n_pi = length(parameter_discretization$pi), 
      tag_list = raw_sattags,
      depth_threshold = deep_dive_depth
    )
  ),
  
  
  tar_target(
    name = nim_fit,
    command = fit(
      nim_pkg = nim_pkg, 
      nsamples = 1e4, 
      nthin = 10, 
      max_batch_iter = 1e3, 
      max_batch_time = 30 * 60, 
      max_batch_mem = 1024^3/2, 
      sample_dir = mcmc_sample_dir
    ),
    deployment = 'worker'
  )
  
  
)
