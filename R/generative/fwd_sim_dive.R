fwd_sim_dive = function(stage, depth_bins, covariates, n_timepoints, nim_pkg,
                        model, lon, lat, depth_threshold, template_bins, 
                        subject_id) {
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
  #  depth_threshold - depth used to generate prop_recent_deep covariate
  #  template_bins - for building out covariates
  #  subject_id - subject number id, so correct tx. mats can be extracted
  
  # initialize stage history
  stages = c()
  
  ind = length(depth_bins)
  n_sampled = 0
  while(n_sampled < n_timepoints) {
    
    #
    # check termination conditions
    #
    
    if(covariates['depth', ind] > depth_threshold) {
      break
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
    cov_inds = seq(
      to = ind + 1,
      from = max(1, ind + 1 - 3600/nim_pkg$consts$tstep)
    )
    covariates_modeled = covariate_tx(
      covariates = covariates[, cov_inds], 
      control = list(
        deep_depth = depth_threshold,
        window_len = 3600,  # TODO: make this an argument somewhere!!!
        obs_freq = nim_pkg$consts$tstep,
        spline_degree = 3,   # TODO: make this an argument somewhere!!!
        lon = lon,
        lat = lat
      )
    )
    
    
    
    #
    # sample stage
    #
    
    # map the new covariate
    covariateId = which(
      apply(nim_pkg$consts$covariates_unique, 2, function(x) {
        all(x == covariates_modeled[, ncol(covariates_modeled)])
      })
    )
    
    # sample next stage
    stage = sample(
      x = 1:nim_pkg$consts$n_stages, 
      size = 1,
      p = model$stage_tx_mat[subject_id, stage, , covariateId]
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
