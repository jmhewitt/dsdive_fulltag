fwd_sim_dive = function(stage, depth_bins, covariates, n_timepoints, nim_pkg,
                        model, lon, lat, depth_threshold, shallow_threshold,
                        template_bins, subject_id, covariate_tx_control) {
  # forward-sampling of dive 
  #
  # Parameters:
  #  stage - dive stage from which to start simulation
  #  depth_bins - sequence of depth bins used to start simulation
  #  covariates - sequence of covariates used to start simulation
  #  n_timepoints - number of timepoints to simulate
  #  nim_pkg - model constants and auxiliary data
  #  model - nimble model containing transition matrices
  #  lon - longitude around which sampling occurs (for day/night calculations)
  #  lat - latitude around which sampling occurs (for day/night calculations)
  #  depth_threshold - sampling ends after this depth is exceeded for "2nd" time
  #  shallow_threshold - simulated dive must go below depth_threshold twice, 
  #    with a visit above shallow_threshold in between in order for simulation
  #    to end
  #  template_bins - for building out covariates
  #  subject_id - subject number id, so correct tx. mats can be extracted
  
  # initialize termination conditions
  visited_deep = FALSE
  recovered_after_deep = FALSE
  
  covariate_tx_control$lon = lon
  covariate_tx_control$lat = lat
  
  # initialize stage history
  stages = c()
  
  ind = length(depth_bins)
  n_sampled = 0
  while(n_sampled < n_timepoints) {
    
    #
    # check termination conditions
    #
    
    ###########################################################################
    #
    #                       termination conditions
    # 
    # simple state machine model to terminate simulation once a recovery 
    # cycle has been completed.  the cycle requires the simulated depth to 
    # exceed the depth threshold, then visit a shallow, recovery depth, then 
    # exceed the depth threshold once more.
    # 
    ###########################################################################
    
    # event: exceeded deep depth
    if(covariates['depth', ind] > depth_threshold) {
      visited_deep = TRUE
      if(recovered_after_deep) {
        break
      }
    } 
    
    # event: in shallow, recovery depth
    if(covariates['depth', ind] < shallow_threshold) {
      if(visited_deep) {
        recovered_after_deep = TRUE
      }
    }
    
    
    #
    # sample
    #
    
    # sample depth bin at next timepoint
    depth_bins[ind+1] = sample(
      x = 1:nim_pkg$consts$n_bins, 
      size = 1,
      prob = model$depth_tx_mat[stage, depth_bins[ind], ]
    )
      
    #
    # update covariates
    #
    
    # extend covariates matrix
    covariates = cbind(covariates, rep(0, nrow(covariates)))
    
    # basic covariates
    covariates['time', ind+1] = covariates['time', ind] + nim_pkg$consts$tstep
    covariates['depth', ind+1] = template_bins$center[depth_bins[ind+1]]
    # only transform covariates wrt. model specification needed for sampling
    nshift = covariate_tx_control$window_len / covariate_tx_control$obs_freq
    cov_inds = seq(
      to = ind + 1,
      from = max(1, ind + 1 - nshift)
    )
    covariates_modeled = covariate_tx(
      covariates = covariates[, cov_inds], 
      control = covariate_tx_control
    )
    
    #
    # sample stage
    #
    
    # # map the new covariate
    # covariateId = which(
    #   apply(nim_pkg$consts$covariates_unique, 2, function(x) {
    #     all(x == covariates_modeled[, ncol(covariates_modeled)])
    #   })
    # )
    
    # # extract the matching transition matrix, or rebuild (i.e., if it is new)
    # if(length(covariateId) == 1) {
    #   txprobs = model$stage_tx_mat[subject_id, stage, , covariateId]
    # } else {
      txprobs = stageTxMats(
        betas = model$beta_tx[subject_id, , , ],
        covariates = covariates_modeled[, ncol(covariates_modeled), drop = FALSE],
        n_timepoints = 1
      )[stage, , ]
    # }
    
    # sample next stage
    stage = sample(
      x = 1:nim_pkg$consts$n_stages, 
      size = 1,
      p = txprobs
    )
    
    # save stage
    stages = c(stages, stage)
    
    #
    # update counters
    # 
    
    n_sampled = n_sampled + 1
    ind = ind + 1
  }
  
  # package results
  if(n_sampled == 0) {
    list(stages = NULL,
         depth_bins = NULL,
         times = NULL)
  } else {
    sampled_inds = seq(
      to = length(depth_bins),
      by = 1,
      length.out = n_sampled
    )
    list(stages = stages, 
         depth_bins = depth_bins[sampled_inds], 
         times = covariates['time', sampled_inds])
  }
}
