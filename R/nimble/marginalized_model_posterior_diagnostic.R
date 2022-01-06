marginalized_model_posterior_diagnostic_script = tar_target(
  name = marginalized_model_posterior_diagnostic,
  command = {
    
    # load posterior samples
    samples = readRDS(fit_marginalized_model$samples)
    
    # get top-level names and groupings of variables being sampled
    sampling_targets = colnames(samples)
    sampling_target_groups = gsub(
      pattern = '\\[.*\\]',
      replacement = '',
      x = sampling_targets
    )
    sampling_groups = unique(sampling_target_groups)
    
    
    # set burn-in
    burn = 1:(nrow(samples)*.1)
    
    #
    # explore posterior output
    #
    
    browser()
    
    summary(effectiveSize(
      mcmc(samples[-burn, sampling_target_groups == 'beta_tx_mu'])
    ))
    
    summary(effectiveSize(
      mcmc(samples[-burn, sampling_target_groups == 'lambda'])
    ))
    
    plot(mcmc(samples[-burn, sampling_target_groups == 'lambda']))
    
    summary(mcmc(samples[-burn,]))
    
    plot(mcmc(samples[-burn, 'lambda[2]']))
    
    0
  }
)