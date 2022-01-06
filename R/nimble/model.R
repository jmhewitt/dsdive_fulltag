modelCode = nimble::nimbleCode({
  
  # prior distributions for speed classes
  for(j in 1:n_lambda_class) {
    lambda[j] ~ dgamma(shape = .01, rate = .01)
  }
  
  # constraint to assist in identifying speed classes (slow, fast)
  constraint_data ~ dconstraint(lambda[1] <= lambda[2])
  
  if(random_effects == TRUE) {
    
    # hierarchical centering for population-level stage effects
    for(i in 1:n_covariates) {
      for(j in 1:n_stages) { # stage being transitioned from
        beta_tx_mu[i,j,1:(n_stages-1)] ~ dmnorm(
          mean = beta_tx_mu_prior_mean[1:(n_stages-1)],
          cov = beta_tx_mu_prior_cov[1:(n_stages-1), 1:(n_stages-1)]
        )
        beta_tx_prec[i,j,1:(n_stages-1),1:(n_stages-1)] ~ dwish(
          R = beta_tx_prec_prior[1:(n_stages-1),1:(n_stages-1)],
          df = beta_tx_prec_prior_k
        )
      }
    }
    
    # individual-level random effects for stage transition matrix covariates
    for(i in 1:n_covariates) {
      for(j in 1:n_stages) { # stage being transitioned from
        for(k in 1:n_subjects) {
          beta_tx[k,i,j,1:(n_stages-1)] ~ dmnorm(
            mean = beta_tx_mu[i,j,1:(n_stages-1)],
            prec = beta_tx_prec[i,j,1:(n_stages-1), 1:(n_stages-1)]
          )
        }
      }
    }
    
  } else {
    
    # population-level stage effects
    for(i in 1:n_covariates) {
      for(j in 1:n_stages) { # stage being transitioned from
        beta_tx_mu[i,j,1:(n_stages-1)] ~ dmnorm(
          mean = beta_tx_mu_prior_mean[1:(n_stages-1)],
          cov = beta_tx_mu_prior_cov[1:(n_stages-1), 1:(n_stages-1)]
        )
      }
    }
    
    # copy population-level effects to all individuals
    for(i in 1:n_covariates) {
      for(j in 1:n_stages) { # stage being transitioned from
        for(k in 1:n_subjects) {
          beta_tx[k,i,j,1:(n_stages-1)] <- beta_tx_mu[i,j,1:(n_stages-1)]
        }
      }
    }
    
  }
  
  # stage transition matrices
  for(k in 1:n_subjects) {
    for(i in 1:n_covariate_combinations) {
      stage_tx_mat[k, 1:n_stages, 1:n_stages, i] <- stageTxMats(
        betas = beta_tx[k, 1:n_covariates, 1:n_stages, 1:(n_stages-1)],
        covariates = covariates_unique[1:n_covariates, i:i],
        n_timepoints = 1
      )[1:n_stages, 1:n_stages, 1]
    }
  }
  
  # construct population-level depth bin transition matrices for each stage
  for(i in 1:n_stages) {
    depth_tx_mat[i, 1:n_bins, 1:n_bins] <- expocall_gpadm(
      H = buildInfinitesimalGenerator(
        pi = pi[stage_defs[i, 1]],
        lambda = lambda[stage_defs[i, 2]],
        M = n_bins,
        stage = 6,
        widths = widths[1:n_bins]
      ),
      t = tstep,
      nrows = n_bins,
      ncols = n_bins
    )[1:n_bins, 1:n_bins]
  }
  
  # largest sampling unit is a sequence of depth bins
  for(seg_num in 1:n_segments) {
    # likelihood for depth bin observations
    depth_bins[segments[seg_num,1]:segments[seg_num,4]] ~ dstatespace(
      obs_lik_dict = depth_tx_mat[1:n_stages, 1:n_bins, 1:n_bins],
      txmat_dict = stage_tx_mat[
        segments[seg_num, 3], 1:n_stages, 1:n_stages, 1:n_covariate_combinations
      ],
      txmat_seq = covariateId[segments[seg_num,1]:segments[seg_num,4]],
      x0 = init_stages[1:n_stages],
      num_obs_states = n_bins,
      num_latent_states = n_stages,
      nt = segments[seg_num,2]
    )
  }
  
})