marginalized_model_posterior_diagnostic_script = tar_target(
  name = marginalized_model_posterior_diagnostic,
  command = {
    
    # tar_load(fit_marginalized_model)
    
    fit_marginalized_model = list(
      samples = file.path('output', 'mcmc', 'fit_marginalized_model'),
      package = file.path('output', 'mcmc', 'fit_marginalized_model', 
                          'nim_pkg.rds')
    )
    
    #
    # load, label, and merge posterior samples
    #
    
    pkg = readRDS(fit_marginalized_model$package)
    
    mvSample_files = dir(
      path = fit_marginalized_model$samples, 
      pattern = 'mvSamples_[0-9]+', 
      full.names = TRUE
    )
    
    mvSample2_files = dir(
      path = fit_marginalized_model$samples, 
      pattern = 'mvSamples2_[0-9]+', 
      full.names = TRUE
    )
    
    samples = do.call(rbind, lapply(mvSample_files, readRDS))
    samples2 = do.call(rbind, lapply(mvSample2_files, readRDS))
    
    colnames(samples) = readRDS(dir(
      path = fit_marginalized_model$samples, 
      pattern = 'mvSamples_colnames',
      full.names = TRUE
    ))
    
    colnames(samples2) = readRDS(dir(
      path = fit_marginalized_model$samples, 
      pattern = 'mvSamples2_colnames',
      full.names = TRUE
    ))
    
    samples = cbind(samples, samples2)
    rm(samples2)
    
    #
    # remove output for quantities that are not actually sampled
    #
    
    # generate node names for fixed effects that are copied within model
    fixef_specs = expand.grid(
      ind = pkg$consts$fixef_inds,
      k = 1:pkg$consts$n_subjects,
      j = 1:pkg$consts$n_stages,
      jm1 = 1:(pkg$consts$n_stages - 1)
    )
    fixef_specs = fixef_specs[complete.cases(fixef_specs),]
    fixef_nodes = apply(fixef_specs, 1, function(r) {
      paste(
        'beta_tx[', 
        paste(r['k'], r['ind'], r['j'], r['jm1'], sep = ', '), 
        ']', 
        sep = ''
      )
    })
    
    # enumerate all of the nodes that are not properly sampled
    not_sampled = c(
      'pi[2]', # treated as known
      fixef_nodes # fixed effects are copied from population-level parameters
    )
    
    not_sampled_inds = which(colnames(samples) %in% not_sampled)
    
    samples = samples[,-not_sampled_inds]
    
    #
    # diagnostics
    #
    
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
    
    plot(mcmc(samples[-burn, sampling_target_groups == 'lambda']))
    
    plot(mcmc(samples[-burn, sampling_target_groups == 'pi']))
    
    for(g in sampling_groups) {
      message(paste('Effective sample size summary for', g))
      print(summary(effectiveSize(
        mcmc(samples[-burn, sampling_target_groups == g])
      )))
      message()
      message(paste('Acceptance rate summary for', g))
      print(summary(1 - rejectionRate(
        mcmc(samples[-burn, sampling_target_groups == g])
      )))
      message('--------------------------------------------------------')
    }
    
    #
    # compare prior and posterior
    #
    
    ggplot(data.frame(x = samples[-burn, 'lambda[1]']), aes(x=x)) +
      stat_density(geom = 'line') + 
      stat_function(fun = function(x) dgamma(x = x, shape = .01, rate = .01), 
                    lty = 3) + 
      xlab(expression(lambda[1])) + 
      ylab(expression(f(lambda[1]))) + 
      theme_few()
      
    0
  }
)