modelCode = nimble::nimbleCode({

  # introduce tx. matrix entries wrt a "dummy" distribution as this is the more 
  # efficient way to introduce large amounts of reference data into the model
  transition_matrices[1:n_txmat_entries] ~ dflatvec(length = n_txmat_entries)
  
  # speed classes
  for(j in 1:n_lambda_class) {
    # continuous parameter
    lambda[j] ~ dgamma(shape = .01, rate = .01)
    # index of discretized parameter
    lambda_ind[j] <- closestIndex(
      value = lambda[j],
      minval = lambda_discretization[1],
      stepsize = lambda_discretization[2],
      nvals = lambda_discretization[3]
    )
  }
  
  # constraint to assist in identifying speed classes
  constraint_data ~ dconstraint(lambda[1] <= lambda[2] & lambda[2] <= lambda[3])
  
  # stage transition matrix
  for(j in 1:n_stages) {
    txmat_stages[j,1:n_stages] ~ ddirch(alpha = prior_alpha[1:n_stages])
  }
  
  # largest sampling unit is a sequence of depth bins
  for(seg_num in 1:n_segments) {
    
    # likelihood for latent stages
    stages[segments[seg_num,1]:segments[seg_num,4]] ~ dstages(
      txmat = txmat_stages[1:n_stages, 1:n_stages],
      n_timepoints = segments[seg_num,2]
    )
    
    # likelihood for depth bin observations
    depths[segments[seg_num,1]:segments[seg_num,4]] ~ dbins(
      # latent stages and definition
      stages = stages[segments[seg_num,1]:segments[seg_num,4]],
      stage_defs = stage_defs[1:n_stages, 1:2],
      # (variable) discretized parameters
      lambda_inds = lambda_ind[1:n_lambda_class],
      # dimension information
      n_timepoints = segments[seg_num,2],
      n_bins = n_bins,
      # transition matrices and discretization
      tmats = transition_matrices[1:n_txmat_entries],
      n_pi = n_pi, 
      n_lambda = lambda_discretization[3]
    )
    
  }
  
})