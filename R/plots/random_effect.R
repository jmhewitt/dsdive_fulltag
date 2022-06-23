random_effect_plot_script = tar_target(
  name = random_effect_plot,
  command = {
    
    #
    # load, label, and merge posterior samples
    #
    
    # tar_load(fit_marginalized_model)
    
    fit_marginalized_model = list(
      samples = file.path('output', 'mcmc', 'fit_marginalized_model'),
      package = file.path('output', 'mcmc', 'fit_marginalized_model', 
                          'nim_pkg.rds')
    )
    
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
    
    pkg = readRDS(fit_marginalized_model$package)
    
    # set burn-in
    burn = 1:(nrow(samples)*.5)
    
    
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
    # summarize random effect posteriors
    #
    
    # generate node names for random effects
    ranef_specs = expand.grid(
      ind = pkg$consts$ranef_inds,
      k = 1:pkg$consts$n_subjects,
      j = 1:pkg$consts$n_stages,
      jm1 = 1:(pkg$consts$n_stages - 1)
    )
    ranef_specs = ranef_specs[complete.cases(ranef_specs),]
    ranef_specs$node = apply(ranef_specs, 1, function(r) {
      paste(
        'beta_tx[', 
        paste(r['k'], r['ind'], r['j'], r['jm1'], sep = ', '), 
        ']', 
        sep = ''
      )
    })
    
    # compute and format posterior summaries for plotting
    df = do.call(rbind, lapply(1:nrow(ranef_specs), function(i) {
      # extract row
      r = ranef_specs[i,]
      # construct hierarchical centering node
      centering_node = paste(
        'beta_tx_mu[',
        paste(r$ind, r$j, r$jm1, sep =', '), 
        ']',
        sep = ''
      )
      # random effect residual wrt. hierarchical centering
      s = samples[-burn, r$node] - samples[-burn, centering_node]
      # package summary
      data.frame(
        ranef_specs[i,],
        mean = mean(s),
        sd = sd(s),
        HPDinterval(mcmc(s)),
        subject = pkg$consts$subject_id_labels[r$k],
        jlab = paste('From: ', rownames(pkg$consts$stage_defs)[r$j], sep =''),
        jm1lab = paste(
          'Offset to log(',
          rownames(pkg$consts$stage_defs)[r$jm1],
          ' / ',
          tail(rownames(pkg$consts$stage_defs),1),
          ')',
          sep = ''
        ),
        label = paste(
          'From: ', 
          rownames(pkg$consts$stage_defs)[r$j],
          ',\nOffset to log(',
          rownames(pkg$consts$stage_defs)[r$jm1],
          ' / ',
          tail(rownames(pkg$consts$stage_defs),1),
          ')',
          sep = ''
        )
      ) %>% mutate(
        signif = !((lower <= 0) & (0 <= upper))
      )
    }))
    
    pl = ggplot(df, aes(x = subject, y = mean, ymin = lower, ymax = upper)) + 
      # non-significant line 
      geom_hline(yintercept = 0, lty = 3) +
      # posterior summaries
      geom_pointrange() + 
      # separate plots for each random effect
      facet_grid(jm1lab~jlab, scales = 'free_x') + 
      # formatting
      ylab('Random effect') + 
      xlab('') +
      coord_flip() + 
      theme_few()
  
    #
    # save plot
    #
    
    f = file.path('output', 'figures')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(pl, filename = file.path(f, paste(tar_name(), '.pdf', sep ='')), 
           dpi = 'print', width = 16, height = 16)
    
  }
)