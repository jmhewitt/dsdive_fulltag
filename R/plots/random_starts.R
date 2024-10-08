random_starts_plot_script = tar_target(
  name = random_starts_plot,
  command = {
    
    #
    # enumerate random starts outputs
    #
    
    outdirs = dir(
      path = file.path('output', 'mcmc', 'fixed_init_beta'), 
      pattern = 'fit_marginalized_model_', 
      full.names = TRUE
    )
    
    # compute posterior summaries for each chain
    df = do.call(rbind, lapply(outdirs, function(d) {
      message(d)
      #
      # load, label, and merge posterior samples
      #
      
      mvSample_files = dir(
        path = d, 
        pattern = 'mvSamples_[0-9]+', 
        full.names = TRUE
      )
      
      mvSample2_files = dir(
        path = d, 
        pattern = 'mvSamples2_[0-9]+', 
        full.names = TRUE
      )
      
      samples = do.call(rbind, lapply(mvSample_files, readRDS))
      samples2 = do.call(rbind, lapply(mvSample2_files, readRDS))
      
      colnames(samples) = readRDS(dir(
        path = d, 
        pattern = 'mvSamples_colnames',
        full.names = TRUE
      ))
      
      colnames(samples2) = readRDS(dir(
        path = d, 
        pattern = 'mvSamples2_colnames',
        full.names = TRUE
      ))
      
      samples = cbind(samples, samples2)
      rm(samples2)
      
      pkg = readRDS(file.path(d, 'nim_pkg.rds'))
      
      # set burn-in
      # burn = 1:(nrow(samples)*.5)
      burn = 1:7e3
      
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
      
      if(nrow(samples) < 1e4) {
        return(NULL)
      }
      
      data.frame(
        rep = as.numeric(str_extract(string = basename(d), pattern = '[0-9]+')),
        param = colnames(samples),
        mean = colMeans(samples[-burn,], na.rm = TRUE),
        t(apply(samples[-burn,], 2, function(x) {
          HPDinterval(mcmc(x))
        }))
      )
    }))
    
    colnames(df)[4:5] = c('lower', 'upper')
    
    
    #
    # make plots for the main parameters sampled
    #
    
    # reminder of what parameters were, indeed, sampled
    sampling_groups = unique(gsub(
      pattern = '\\[.*\\]', replacement = '', x = unique(df$param)
    ))
      
    # posteriors for covariates that influence depth bin transitions (or at 
    # least just the intercept)
    pl = ggplot(
      df %>% filter(grepl(pattern = 'beta_tx_mu\\[1', x = param)),
      aes(x = rep, y = mean, ymin = lower, ymax = upper)
    ) + 
      geom_pointrange() + 
      facet_wrap(~param, scales = 'free_y') + 
      theme_few()
    
    # posteriors for speed parameters
    pl2 = ggplot(
      df %>% filter(grepl(pattern = 'lambda', x = param)),
      aes(x = rep, y = mean, ymin = lower, ymax = upper)
    ) + 
      geom_pointrange() + 
      facet_wrap(~param, scales = 'free_y') + 
      theme_few()
    
    # posteriors for directional bias parameters
    pl3 = ggplot(
      df %>% filter(grepl(pattern = 'pi', x = param)),
      aes(x = rep, y = mean, ymin = lower, ymax = upper)
    ) + 
      geom_pointrange() + 
      facet_wrap(~param, scales = 'free_y') + 
      theme_few()
    
    #
    # save plots
    #
    
    f = file.path('output', 'figures', 'parameter_sensitivity')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(pl, 
           filename = file.path(f, paste(tar_name(), '_beta_tx_mu.pdf', 
                                         sep = '')),
           width = 16, height = 16)
    
    ggsave(pl2, 
           filename = file.path(f, paste(tar_name(), '_lambda.pdf', 
                                         sep = '')))
    
    ggsave(pl3, 
           filename = file.path(f, paste(tar_name(), '_pi.pdf', 
                                         sep = '')))
    
  }
)