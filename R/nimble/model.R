modelCode = nimbleCode({

  # introduce tx. matrix entries wrt a "dummy" distribution as this is the more 
  # efficient way to introduce large amounts of reference data into the model
  transition_matrices[1:n_txmat_entries] ~ dflatvec()
  
  for(i in 1:n_subjects) {
    betas_tx[1:n_covariates, 1:11, i] ~ dflatmat()
  }
  
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

  for(seg_num in 1:n_segments) {
    
    # introduce covariates wrt a "dummy" distribution as this is the more 
    # efficient way to introduce large numbers of covariates into the model
    covariates[
      1:n_covariates,  segments[seg_num,1]:segments[seg_num,4]
    ] ~ dflatmat()
    
    stages[segments[seg_num,1]:segments[seg_num,4]] ~ dflatvec()
    
    depths[segments[seg_num,1]:segments[seg_num,4]] ~ dbins(
      # latent stages
      stages = stages[segments[seg_num,1]:segments[seg_num,4]], 
      # dimension information
      n_timepoints = segments[seg_num,2], 
      stage_map = movement_types[1:n_stages], n_bins = n_bins, 
      # depth bin transition covariates for the subject associated w/segment
      alpha = alpha[1:n_covariates, 1:n_stages, segments[seg_num, 3]],
      beta = beta[1:n_covariates, 1:n_stages, segments[seg_num, 3]],
      covariates = covariates[
        1:n_covariates,  segments[seg_num,1]:segments[seg_num,4]
      ],
      # transition matrices and discretization
      tmats = transition_matrices[1:n_txmat_entries],
      pi_discretization = pi_discretization[1:n_txmat_types, 1:3],
      n_pi = n_pi[1:n_txmat_types], n_lambda = n_lambda[1:n_txmat_types],
      lambda_discretization = lambda_discretization[1:n_txmat_types, 1:3]
    )
    
  }
  
})