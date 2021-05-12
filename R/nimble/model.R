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
  
  # stage transition matrix effects
  for(i in 1:n_covariates) {
    for(j in 1:n_stages) { # stage being transitioned from
      beta_tx[i,j,1:(n_stages-1)] ~ dmnorm(mean = beta_tx_mu[1:(n_stages-1)],
                                           cov = beta_tx_cov[1:(n_stages-1),
                                                             1:(n_stages-1)])
    }
  }
  
  # largest sampling unit is a sequence of depth bins
  for(seg_num in 1:n_segments) {
    
    # introduce covariates wrt a "dummy" distribution as this is the more
    # efficient way to introduce large numbers of covariates into the model
    covariates[
      1:n_covariates,  segments[seg_num,1]:segments[seg_num,4]
    ] ~ dflatmat(nrow = n_covariates, 
                 ncol = segments[seg_num,4] - segments[seg_num,1] + 1)
    
    # likelihood for latent stages
    stages[segments[seg_num,1]:segments[seg_num,4]] ~ dstages(
      beta_tx = beta_tx[1:n_covariates, 1:n_stages, 1:(n_stages-1)],
      covariates = covariates[1:n_covariates,  
                              segments[seg_num,1]:segments[seg_num,4]],
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