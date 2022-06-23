random_starts_plot_script = tar_target(
  name = random_starts_plot,
  command = {
    
    browser()
    
    #
    # enumerate random starts outputs
    #
    
    outdirs = dir(
      path = file.path('output', 'mcmc'), 
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
    
    
    pl = ggplot(
      df %>% filter(grepl(pattern = 'beta_tx_mu\\[1', x = param)),
      aes(x = rep, y = mean, ymin = lower, ymax = upper)
    ) + 
      geom_pointrange() + 
      facet_wrap(~param, scales = 'free_y') + 
      theme_few()
  
    pl  
  }
)