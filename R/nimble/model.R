modelCode = nimbleCode({

  # introduce tx. matrix entries wrt a "dummy" distribution as this is the more 
  # efficient way to introduce large amounts of reference data into the model
  transition_matrices[1:n_txmat_entries] ~ dflatvec(length = n_txmat_entries)
  
  # TODO: give a proper prior and population-level model to stage transitions!
  # TODO: add covariates s.t. we only see certain types of transitions at surface
  for(i in 1:n_subjects) {
    for(j in 1:n_covariates) {
      for(k in 1:n_stage_txs) {
        betas_tx[j, k, i] ~ dnorm(mean = 0, sd = 1e2)
      }
    }
  }
  
  # TODO: soft-code hyper priors to allow sensitivity studies, etc.
  # depth bin transition model random effects and population-level parameters
  for(i in 1:n_covariates) {
    for(j in 1:n_stages) {
      # priors for population-level coefficients and variability
      alpha_mu[i,j] ~ dnorm(mean = 0, sd = 1e2)
      beta_mu[i,j] ~ dnorm(mean = 0, sd = 1e2)
      alpha_var[i,j] ~ dinvgamma(shape = 2, rate = 1)
      beta_var[i,j] ~ dinvgamma(shape = 2, rate = 1)
      # hierarchical layer for individual-level coefficients
      for(k in 1:n_subjects) {
        alpha[i,j,k] ~ dnorm(mean = alpha_mu[i,j], var = alpha_var[i,j])
        beta[i,j,k] ~ dnorm(mean = beta_mu[i,j], var = beta_var[i,j])
      }
    }
  }

  # largest sampling unit is a sequence of depth bins
  for(seg_num in 1:n_segments) {
    
    # # introduce covariates wrt a "dummy" distribution as this is the more 
    # # efficient way to introduce large numbers of covariates into the model
    covariates[
      1:n_covariates,  segments[seg_num,1]:segments[seg_num,4]
    ] ~ dflatmat(nrow = n_covariates, 
                 ncol = segments[seg_num,4] - segments[seg_num,1] + 1)
    
    # TODO: give a proper likelihood to stages!
    # stages[segments[seg_num,1]:segments[seg_num,4]] ~ dflatvec()
    
    # likelihood for depth bin observations
    depths[segments[seg_num,1]:segments[seg_num,4]] ~ dbins(
      # latent stages
      stages = stages[segments[seg_num,1]:segments[seg_num,4]],
      # dimension information
      n_bins = n_bins,
      n_timepoints = segments[seg_num,2],
      stage_map = movement_types[1:n_stages],
      # depth bin transition covariates for the subject associated w/segment
      covariates = covariates[
        1:n_covariates,  segments[seg_num,1]:segments[seg_num,4]
      ],
      beta = beta[1:n_covariates, 1:n_stages, segments[seg_num, 3]],
      alpha = alpha[1:n_covariates, 1:n_stages, segments[seg_num, 3]],
      # transition matrices and discretization
      n_pi = n_pi[1:n_txmat_types],
      n_lambda = n_lambda[1:n_txmat_types],
      tmats = transition_matrices[1:n_txmat_entries],
      pi_discretization = pi_discretization[1:n_txmat_types, 1:3],
      lambda_discretization = lambda_discretization[1:n_txmat_types, 1:3]
    )
    
  }
  
})