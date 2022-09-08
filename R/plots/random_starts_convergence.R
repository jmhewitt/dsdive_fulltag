random_starts_convergence_script = tar_target(
  name = random_starts_convergence,
  command = {
    
    #
    # enumerate random starts outputs
    #
    
    outdirs = dir(
      path = file.path('output', 'mcmc', 'fixed_init_beta'), 
      pattern = 'fit_marginalized_model_', 
      full.names = TRUE
    )
    
    # aggregate parameter samples from each chain
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
      
      if(length(mvSample_files) == 0) {
        return(NULL)
      }
      
      samples = do.call(rbind, lapply(mvSample_files, readRDS))
      
      colnames(samples) = readRDS(dir(
        path = d, 
        pattern = 'mvSamples_colnames',
        full.names = TRUE
      ))
      
      # define sampling targets
      tgt = c(
        paste('pi[',c(1,3),']', sep =''),
        paste('lambda[',1:2,']', sep ='')
      )
      
      samples = samples[,tgt]
      
      if(nrow(samples) < 1e4) {
        return(NULL)
      }
      
      data.frame(
        samples,
        ind = 1:nrow(samples),
        rep = d
      )
    }))
    
    # format parameter names for plotting
    df_longer = df %>% 
      pivot_longer(cols = colnames(df)[1:4]) %>% 
      mutate(
        name = stringr::str_replace(
          string = name, 
          pattern = '\\.([0-9]+)\\.', 
          replacement = '[\\1]'
        )
      )
    
    # plot beginnings of chains
    pl_begin = ggplot(df_longer %>% filter(1 <= ind, ind <= 1e3), 
                      aes(x = ind, y = value, group = rep)) + 
      geom_line(alpha = .2) + 
      facet_wrap(~name, scales = 'free_y', labeller = label_parsed) + 
      guides(col = 'none') + 
      xlab('Iteration') + 
      ylab(label = 'Value') + 
      theme_few()
    
    # plot ends of chains
    pl_end = ggplot(df_longer %>% filter(8.5e3 <= ind, ind <= 10e3), 
                      aes(x = ind, y = value, group = rep)) + 
      geom_line(alpha = .2) + 
      facet_wrap(~name, scales = 'free_y', labeller = label_parsed) + 
      guides(col = 'none') + 
      xlab('Iteration') + 
      ylab(label = 'Value') + 
      theme_few()
    
    # create output directory
    f = file.path('output', 'figures', 'convergence')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(
      pl_begin, 
      filename = file.path(f, paste(tar_name(), '_multiple_chains_start.pdf', 
                                    sep = '')), 
      dpi = 'print', width = 8, height = 6
    )
    
    ggsave(
      pl_end, 
      filename = file.path(f, paste(tar_name(), '_multiple_chains_end.pdf', 
                                    sep = '')), 
      dpi = 'print', width = 8, height = 6
    )
    
    
    # densities
    pl_densities = ggplot(df_longer %>% filter(8.5e3 <= ind, ind <= 10e3), 
                aes(x = value)) + 
      stat_density(geom = 'line') + 
      facet_wrap(~name, labeller = label_parsed, scales = 'free', 
                 strip.position = 'bottom') + 
      xlab('Parameter') + 
      ylab('Posterior density') + 
      theme_few() + 
      theme(strip.placement = 'outside', 
            axis.title.x = element_blank())
    
    ggsave(
      pl_densities, 
      filename = file.path(f, paste(tar_name(), '_combined_posteriors.pdf', 
                                    sep = '')), 
      dpi = 'print', width = 8, height = 6
    )
    
    # density summaries
    sink(file.path(f, paste(tar_name(), 
                            '_combined_posterior_summaries.txt', sep = '')))
    df_filtered = df %>% filter(8.5e3 <= ind, ind <= 10e3)
    data.frame(
      mean = round(colMeans(mcmc(df_filtered[,1:4])), 2),
      sd = round(apply(df_filtered[,1:4], 2, sd), 3),
      round(HPDinterval(mcmc(df_filtered[,1:4])), 2)
    )
    sink()
    
    0
  }
)