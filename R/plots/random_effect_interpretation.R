random_effect_interpretation_plot_script = tar_target(
  name = random_effect_interpretation_plot,
  command = {
    
    # identify posterior predictive samples
    f = dir(
      path = file.path(
        'output', 'parameter_interpretation', 'fixed_init_beta', 
        'random_effects', 'samples'
      ), 
      full.names = TRUE
    )
    
    # get subject id labels
    pkg = readRDS(
      file.path('output', 'mcmc', 'fit_marginalized_model_0', 'nim_pkg.rds')
    )
    
    # load posterior samples
    samples = do.call(rbind, lapply(f, function(f) {
      r = readRDS(f)[[1]]
      cbind(
        data.frame(r$config),
        # seconds to minutes
        data.frame(t(r$baseline_deep_pred_samples / 60))
      )
    }))
    
    # compute posterior summaries
    summaries = do.call(rbind, lapply(f, function(f) {
      r = readRDS(f)[[1]]
      # seconds to minutes
      r$baseline_deep_pred_samples = r$baseline_deep_pred_samples / 60
      # extract samples
      h1 = r$baseline_deep_pred_samples['first_deep',]
      h2 = r$baseline_deep_pred_samples['second_deep',] - h1
      # summaries
      cbind(
        data.frame(r$config),
        data.frame(
          parameter = c('first_deep','recovery'),
          mean = c(mean(h1), mean(h2)),
          sd = c(sd(h1), sd(h2))
        ),
        rbind(
          HPDinterval(mcmc(h1)),
          HPDinterval(mcmc(h2))
        )
      )
    }))
    
    # compute posterior summaries
    densities = do.call(rbind, lapply(f, function(f) {
      r = readRDS(f)[[1]]
      # seconds to minutes
      r$baseline_deep_pred_samples = r$baseline_deep_pred_samples / 60
      # extract samples
      h1 = r$baseline_deep_pred_samples['first_deep',]
      h2 = r$baseline_deep_pred_samples['second_deep',] - h1
      # summaries
      d1 = density(h1)
      d2 = density(h2)
      cbind(
        data.frame(r$config),
        rbind(
          data.frame(
            parameter = 'first_deep',
            value = d1$x,
            density = d1$y
          ),
          data.frame(
            parameter = 'recovery',
            value = d2$x,
            density = d2$y
          )
        )
      )
    }))
    
    # label the different trajectories used to initialize simulations
    summaries$prototype = factor(
      x = summaries$prototype, 
      levels = rev(1:3), 
      labels = rev(c(
        'Shallow', 'Recovery', 'Deep'
      ))
    )
    
    # label the different trajectories used to initialize simulations
    densities$prototype = factor(
      x = densities$prototype, 
      levels = rev(1:3), 
      labels = rev(c(
        'Shallow', 'Recovery', 'Deep'
      ))
    )
    
    # relabel the lighting conditions
    summaries$time = factor(
      x = summaries$time,
      levels = c("daytime", "night_dark", "night_moonlit")[c(1,3,2)],
      # labels = c('Daytime', 'Dark\nNight', 'Moonlit\nNight')[c(1,3,2)]
      labels = c('Daytime', 'Dark Night', 'Moonlit Night')[c(1,3,2)]
    )
    
    # relabel the lighting conditions
    densities$time = factor(
      x = densities$time,
      levels = c("daytime", "night_dark", "night_moonlit")[c(1,3,2)],
      # labels = c('Daytime', 'Dark\nNight', 'Moonlit\nNight')[c(1,3,2)]
      labels = c('Daytime', 'Dark Night', 'Moonlit Night')[c(1,3,2)]
    )
    
    # convert inputs like 'night_dark' to 'Night Dark'
    titlecase = function(labels) {
      .simpleCap = function(x) {
        sapply(x, function(x) {
          s = strsplit(x, " ")[[1]]
          paste(toupper(substring(s, 1, 1)), substring(s, 2),
                sep = "", collapse = " ")
        })
      }
      lapply(labels, function(lab) {
        
        .simpleCap(gsub('_', ' ', as.character(lab)))
      })
    }
    
    # visualize random effects
    pl = ggplot(summaries, aes(x = mean, y = sd)) + 
      geom_point() + 
      facet_grid(parameter~prototype, scales = 'free') + 
      theme_few()
    
    pl = ggplot(summaries, aes(x = subject_id, y = mean, ymin = lower, ymax = upper)) + 
      geom_pointrange() + 
      facet_grid(parameter~prototype, scales = 'free') + 
      theme_few()
    
    pl = ggplot(samples, aes(x = factor(subject_id), y = first_deep)) + 
      geom_boxplot() + 
      facet_wrap(~prototype, scales = 'free') + 
      theme_few()
    
    pl = ggplot(densities %>% filter(parameter == 'first_deep'), 
                aes(x = value, y = density, group = subject_id)) + 
      geom_line(alpha = .1) + 
      facet_wrap(~prototype, scales = 'free') + 
      xlim(0, 300) +
      xlab(expression(H[1](1))) + 
      ylab('Post. Pred. density') + 
      theme_few()
    
    #
    # save figures
    #
    
    f = file.path('output', 'figures', 'fixed_init_beta')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(pl, 
           filename = file.path(f, paste(tar_name(), '.pdf', sep = '')),
           width = 8, height = 3, dpi = 'print')
    
    #
    # save raw summaries
    #
    
    summaries$subject_id = pkg$consts$subject_id_labels[summaries$subject_id]
    
    saveRDS(
      summaries %>% filter(parameter == 'first_deep'), 
      file = file.path(f, paste(tar_name(), '_data.rds', sep=''))
    )
      
    summaries %>% 
      filter(parameter == 'first_deep') %>% 
      select(sd) %>% unlist() %>% range() %>% round()
    
  }, 
  deployment = 'worker', 
  memory = 'transient',
  storage = 'worker'
)