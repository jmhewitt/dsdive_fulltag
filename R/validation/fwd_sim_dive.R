fwd_sim_dive = function(stages, depths, covariates, n_timepoints, nim_pkg,
                        alpha, beta, betas_tx, template_bins, times, timestep, 
                        lon, lat) {
  # forward-sampling of dive 
  #
  # Parameters:
  #  stages - dive stage from which to start simulation
  #  depths - depth bin from which to start simulation
  #  covariates - covariates with which to start simulation
  #  n_timepoints - number of timepoints to simulate, or Inf if dive should be 
  #    simulated until a surface bin is reached
  #  nim_pkg - model constants
  #  alpha - descent preference parameters
  #  beta - dive speed parameters
  #  betas_tx - stage transition parameters
  #  times - time at which sampling starts
  #  timestep - time between observations/samples
  #  lon - longitude around which sampling occurs (for day/night calculations)
  #  lat - latitude around which sampling occurs (for day/night calculations)
  
  ind = 1
  while(ind < n_timepoints) {
    
    #
    # sample depth bin at next timepoint
    #
    
    # extract movement type associated with dive stage
    s <- stages[ind]
    mtype <- nim_pkg$consts$movement_types[s]
    # compute pi and lambdas
    pi_val <- ilogit(inprod(alpha[,s], covariates[,ind]))
    lambda_val <- exp(inprod(beta[,s], covariates[,ind]))
    # discretize pi and lambdas
    pi_ind  <- closestIndex(
      value = pi_val,
      minval = nim_pkg$consts$pi_discretization[mtype, 1],
      stepsize = nim_pkg$consts$pi_discretization[mtype, 2],
      nvals = nim_pkg$consts$pi_discretization[mtype, 3]
    )
    lambda_ind  <- closestIndex(
      value = lambda_val,
      minval = nim_pkg$consts$lambda_discretization[mtype, 1],
      stepsize = nim_pkg$consts$lambda_discretization[mtype, 2],
      nvals = nim_pkg$consts$lambda_discretization[mtype, 3]
    )
    # relevant row from discretized transition matrix
    p = lookupTmat(
      movement_type = mtype, pi_ind = pi_ind, lambda_ind = lambda_ind,
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
    covariates['depth',ind+1] = template_bins$center[depths[ind+1]]
    covariates['deep_depth',ind+1] = covariates['depth',ind+1] > 800
    covariates['shallow_depth',ind+1] = covariates['depth',ind+1] <= 800
    covariates['non_surface_bin',ind+1] = depths[ind+1] > 1
    covariates['surface_bin',ind+1] = depths[ind+1] == 1
    
    # event-based covariates
    covariates['time_since_surface',ind+1] =  ifelse(
      covariates['surface_bin',ind+1] == TRUE,
      0, 1 + covariates['time_since_surface',ind]
    )
    covariates['all_shallow_depths_since_surface',ind+1] =  ifelse(
      covariates['all_shallow_depths_since_surface',ind] == FALSE,
      FALSE, ifelse(covariates['shallow_depth',ind+1] == TRUE, TRUE, FALSE)
    )
    
    # celestial covariates
    times[ind + 1] = times[ind] + timestep
    covariates['daytime',ind+1] = daytime(date = times[ind+1], lat = lat, 
                                          lon = lon)
    covariates['moonlit',ind+1] = moonlit(date = times[ind+1], lat = lat, 
                                          lon = lon)

    
    #
    # sample stage
    #
    
    # relevant row from transition matrix
    p = stageTxVec(stageFrom = s, betas = betas_tx, 
                   covariates = covariates[,ind+1], 
                   surface_bin = covariates['surface_bin',ind+1],
                   log = FALSE)
    
    # next stage
    stages[ind+1] = sample(1:nim_pkg$consts$n_stages, size = 1, prob = p)
    
    
    #
    # check termination conditions and loop
    #
    
    if(is.infinite(n_timepoints)) {
      # stop sampling at surface
      if(covariates['surface_bin',ind+1] == TRUE) {
        break
      }
    }
    
    # increment counter
    ind = ind + 1
  }
  
  # package results
  list(stages = stages, depths = depths, times = times)
}