fwd_sim_to_depth = function(stages, depths, covariates, n_max, nim_pkg,
                           lambda, betas_tx, template_bins, times, timestep, 
                           lon, lat, depth_threshold, deeper = TRUE) {
  # forward-sampling of dive 
  #
  # Parameters:
  #  stages - dive stage from which to start simulation
  #  depths - depth bin from which to start simulation
  #  covariates - covariates with which to start simulation
  #  n_max - maximum number of timepoints to simulate
  #  nim_pkg - model constants
  #  lambdas - speed parameters
  #  betas_tx - stage transition parameters
  #  times - time at which sampling starts
  #  timestep - time between observations/samples
  #  lon - longitude around which sampling occurs (for day/night calculations)
  #  lat - latitude around which sampling occurs (for day/night calculations)
  #  depth_threshold - depth used to generate prop_recent_deep covariate
  #  deeper - TRUE if depth_threshold should be exceeded
  
  # discretize speed parameters
  lambda_inds = sapply(lambda, function(par) {
    closestIndex(
      value = par,
      minval = nim_pkg$consts$lambda_discretization[1],
      stepsize = nim_pkg$consts$lambda_discretization[2],
      nvals = nim_pkg$consts$lambda_discretization[3]
    )
  })
  
  # save this so we can filter the return output
  start_ind = length(depths)
  
  ind = start_ind
  n_sampled = 0
  while(n_sampled < n_max) {
    
    #
    # sample depth bin at next timepoint
    #
    
    # extract pi/lambda values
    pi_ind = nim_pkg$consts$stage_defs[stages[ind], 1]
    lambda_ind = lambda_inds[nim_pkg$consts$stage_defs[stages[ind], 2]]
    # relevant row from discretized transition matrix
    p = lookupTmat(
      pi_ind = pi_ind, lambda_ind = lambda_ind,
      n_pi = nim_pkg$consts$n_pi, n_lambda = nim_pkg$consts$n_lambda,
      n_bins = nim_pkg$consts$n_bins, tmats = nim_pkg$data$transition_matrices
    )[depths[ind],]
    # next depth
    depths[ind+1] = sample(1:nim_pkg$consts$n_bins, size = 1, prob = p)
      
    #
    # update covariates
    #
    
    # extend covariates matrix
    covariates = cbind(covariates, rep(0, nrow(covariates)))
    
    # basic covariates
    covariates['intercept',ind+1] = 1
    
    # celestial covariates
    times[ind + 1] = times[ind] + timestep
    covariates['daytime',ind+1] = daytime(date = times[ind+1], lat = lat, 
                                          lon = lon)
    covariates['moonlit',ind+1] = moonlit(date = times[ind+1], lat = lat, 
                                          lon = lon)
    
    #
    # proportion of recent observations at depth, as a covariate
    #
      
    # window at which recent observations begins
    window_start = times[ind] - duration(1, units = 'hours')
    # data indices of recent observations
    past_inds = 1:(ind-1)
    window_inds = past_inds[window_start <= times[past_inds]]
    # compute proportion of recent observations spent below a depth
    covariates['prop_recent_deep',ind+1] =  sum(
      template_bins$center[depths[window_inds]] >= depth_threshold
    ) / length(window_inds)
    covariates['prop_recent_deep3',ind+1] =  
      (covariates['prop_recent_deep',ind+1] - .5)^3
    
    
    #
    # sample stage
    #
    
    
    # relevant row from stage transition matrix
    p = stageTxVec(stageFrom = stages[ind], betas = betas_tx, 
                   covariates = covariates[,ind+1], log = FALSE)
    
    # next stage
    stages[ind+1] = sample(1:nim_pkg$consts$n_stages, size = 1, prob = p)
    
    
    #
    # check termination conditions and loop
    #
    
    
    if(deeper) {
      # stop sampling when a deep depth is reached
      if(template_bins$center[depths[ind + 1]] >= depth_threshold) {
        break
      }
    } else {
      # stop sampling when a shallow depth is reached
      if(template_bins$center[depths[ind + 1]] <= depth_threshold) {
        break
      }
    }
    
    # increment counter
    n_sampled = n_sampled + 1
    ind = ind + 1
  }
  
  len = length(stages)
  
  # package results
  list(stages = stages[start_ind:len], 
       depths = depths[start_ind:len], 
       times = times[start_ind:len])
}