
#
# dive model
#

modelCode_stageLearning = nimbleCode({
  
  #
  # priors
  #
  
  # tag-specific stage duration parameters
  for(i in 1:N_tags) {
    
    # loop over deep dive stages
    for(j in 1:2) {
      # marginal means
      xi_prior_means[i, j] ~ dnorm(mean = G_prior_mean[j], sd = G_prior_sd[j])
      # marginal variances
      xi_prior_covs[i, j, j] ~ dgamma(shape = G_prior_shape[j],
                                      rate = G_prior_rate[j])
    }
    
    # deep dive: bivariate correlation
    xi_prior_cor_scaled[i] ~ dbeta(shape1 = G_prior_cor[1], 
                                   shape2 = G_prior_cor[2])
    xi_prior_cor[i] <- xi_prior_cor_scaled[i] * 2 - 1
    # induced off-diagonal covariance
    xi_prior_covs[i, 1, 2] <- xi_prior_cor[i] * 
      sqrt(xi_prior_covs[i, 1, 1]) * 
      sqrt(xi_prior_covs[i, 2, 2])
    xi_prior_covs[i, 2, 1] <- xi_prior_covs[i, 1, 2]
    
    # shallow dives: marginal mean
    xi_shallow_prior[i, 1] ~ dnorm(mean = G_shallow_prior_mean[1], 
                                   sd = G_shallow_prior_mean[2])
    # marginal variance
    xi_shallow_prior[i, 2] ~ dgamma(shape = G_shallow_prior_sd[1],
                                    rate = G_shallow_prior_sd[2])
  }
  
  # population level directional preferences
  for(i in c(1,3,4,5)) {
    logit_pi[i] ~ dlogitBeta(shape1 = pi_prior[i, 1], shape2 = pi_prior[i, 2])
    pi[i] <- ilogit(logit_pi[i])
  }
  
  # population level speeds
  for(i in 1:5) {
    for(j in 1:2) {
      beta_lambda[i, j] ~ dnorm(mean = beta_lambda_prior_mean[i, j], 
                                sd = beta_lambda_prior_sd[i, j])
    }
    sigma_lambda[i] ~ dinvgamma(shape = sigma_lambda_priors[i, 1],
                                rate = sigma_lambda_priors[i, 2])
  }
  
  # tag-level movement parameters
  for(tagInd in 1:N_tags) {
    
    # speeds and CTMC generator decompositions
    for(i in 1:5) {
      log_lambda[tagInd, i] ~ dnorm(
        mean = beta_lambda[i, 1] + 
          beta_lambda[i, 2] * tag_covariates[tagInd, 1],
        var = sigma_lambda[i]
      )
      
      # transformed parameters, for sampling
      lambda[tagInd, i] <- exp(log_lambda[tagInd, i])
      
      # generators
      expm_decomp[tagInd, i, 1:(2*N_bins+3), 1:N_bins] <- 
        buildAndDecomposeGenerator(
          pi = pi[i], lambda = lambda[tagInd, i], M = N_bins, stage = i,
          widths = widths[1:N_bins], delta = delta, t = tstep
        )
    }
    
  }
  
  # dive endpoint random effects
  for(i in 1:N_endpoints) {
    endpoints[i] ~ dunif(endpoint_priors[i, 1], endpoint_priors[i, 2])
  }
  
  #
  # likelihood, and additional priors
  #

  # deep dives
  for(diveId in 1:N_dives) {

    # stage duration priors
    log_xi[diveId, 1:2] ~ dmnorm(
      mean = xi_prior_means[dive_relations[diveId, 3], 1:2], 
      cov = xi_prior_covs[dive_relations[diveId, 3], 1:2, 1:2]
    )
    
    # back-transform stage durations
    xi[diveId, 1] <- exp(log_xi[diveId, 1])
    xi[diveId, 2] <- exp(log_xi[diveId, 2])
    
    # dive start and end times
    T[diveId, 1] <- endpoints[dive_relations[diveId, 1]]
    T[diveId, 4] <- endpoints[dive_relations[diveId, 2]]
    
    # stage transition times
    T[diveId, 2] <- T[diveId, 1] + xi[diveId, 1]
    T[diveId, 3] <- T[diveId, 2] + xi[diveId, 2]

    # likelihood
    depths[dive_relations[diveId, 4]:dive_relations[diveId, 5]] ~ ddive(
      times = times[dive_relations[diveId, 4]:dive_relations[diveId, 5]],
      N = dive_relations[diveId, 5] - dive_relations[diveId, 4] + 1,
      expm = expm_decomp[dive_relations[diveId, 3], 1:3, 1:N_bins, 1:N_bins],
      evecs = expm_decomp[dive_relations[diveId, 3], 1:3,
                          (N_bins + 1):(2 * N_bins), 1:N_bins],
      evals = expm_decomp[dive_relations[diveId, 3], 1:3, 2 * N_bins + 1,
                          1:N_bins],
      d = expm_decomp[dive_relations[diveId, 3], 1:3, 2 * N_bins + 2, 1:N_bins],
      dInv = expm_decomp[dive_relations[diveId, 3], 1:3, 2* N_bins + 3,
                         1:N_bins],
      tstep = tstep, M = N_bins, T = T[diveId, 1:4], include = 1
    )
  }
  
  # shallow dives
  for(diveId in 1:N_dives_shallow) {
    
    # stage duration priors
    log_xi_shallow[diveId] ~ dnorm(
      mean = xi_shallow_prior[dive_relations_shallow[diveId, 3], 1], 
      var = xi_shallow_prior[dive_relations_shallow[diveId, 3], 2]
    )
    
    # back-transform stage durations
    xi_shallow[diveId] <- exp(log_xi_shallow[diveId])
    
    # dive start and end times
    T_shallow[diveId, 1] <- endpoints[dive_relations_shallow[diveId, 1]]
    T_shallow[diveId, 3] <- endpoints[dive_relations_shallow[diveId, 2]]
    
    # stage transition times
    T_shallow[diveId, 2] <- T_shallow[diveId, 1] + xi_shallow[diveId]
    
    # likelihood
    depths[
      dive_relations_shallow[diveId, 4]:dive_relations_shallow[diveId, 5]
    ] ~ ddiveShallow(
      times = times[
        dive_relations_shallow[diveId, 4]:dive_relations_shallow[diveId, 5]
      ],
      N = dive_relations_shallow[diveId, 5] - 
        dive_relations_shallow[diveId, 4] + 1,
      expm = expm_decomp[
        dive_relations_shallow[diveId, 3], 4:5, 1:N_bins, 1:N_bins
      ],
      evecs = expm_decomp[dive_relations_shallow[diveId, 3], 4:5,
                          (N_bins + 1):(2 * N_bins), 1:N_bins],
      evals = expm_decomp[dive_relations_shallow[diveId, 3], 4:5, 
                          2 * N_bins + 1, 1:N_bins],
      d = expm_decomp[dive_relations_shallow[diveId, 3], 4:5, 2 * N_bins + 2, 
                      1:N_bins],
      dInv = expm_decomp[dive_relations_shallow[diveId, 3], 4:5, 2* N_bins + 3,
                         1:N_bins],
      tstep = tstep, M = N_bins, T = T_shallow[diveId, 1:3], include = 1
    )
  }
  
})
