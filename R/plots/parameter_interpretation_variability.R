parameter_interpretation_variability_plot_script = tar_target(
  name = multiple_start_reps,
  command = {
    
    parameter_interp_rep = paste(
      'parameter_interpretation_', multiple_start_reps, sep =''
    )
    
    # identify posterior predictive samples
    f = dir(
      path = file.path(
        'output', 'parameter_interpretation', parameter_interp_rep, 'samples'
      ), 
      full.names = TRUE, 
      pattern = 'cfg'
    )
    
    # compute means for all groups from posterior predictive samples
    summaries = do.call(rbind, lapply(f, function(f) {
      samples = readRDS(f)[[1]]
      data.frame(
        samples$config,
        first_deep = samples$baseline_deep_pred_samples['first_deep',]
      )
    }))
    
    # label the different trajectories used to initialize simulations
    summaries$prototype = factor(
      x = summaries$prototype, 
      levels = rev(1:3), 
      labels = rev(c(
        # 'mid_recovery','starting_recovery','finishing_deep'
        # 'Low','Medium','High'
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
    
    # relabel the random walk stage
    summaries$init_stage = gsub(
      pattern = 'slow_free',
      replacement = 'slow_random', 
      x = summaries$init_stage
    )
    
    # convert all seconds to minutes
    summaries$first_deep = summaries$first_deep / 60
    # summaries$second_deep = summaries$second_deep / 60
    # summaries$recovery = summaries$recovery / 60
    
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
    
    # alternate visualization for parameter effects
    plalt = ggplot(summaries, aes(x = first_deep, col = init_stage, 
                               lty = init_stage, pch = init_stage,
                               group = init_stage)) + 
      stat_density(geom = 'line') + 
      scale_linetype_discrete('Initial movement', labels = titlecase) + 
      scale_shape_discrete('Initial movement', labels = titlecase) + 
      scale_color_brewer('Initial movement', labels = titlecase, 
                         type = 'qual', palette = 'Dark2') + 
      # xlab('Recent diving activity') + 
      xlim(0,400) + 
      # scale_x_log10(
      #   minor_breaks = c(
      #     seq(from = 1, to = 10, by = 1),
      #     seq(from = 10, to = 100, by = 10),
      #     seq(from = 100, to = 1000, by = 100)
      #   )
      # ) + 
      ylab('Post. Pred. density') +
      xlab(expression(H[1](1))) + 
      facet_grid(prototype~time, labeller = titlecase) + 
      theme_few() # + 
      # theme(
      #   panel.grid.major.x = element_line(color = 'grey80'),
      #   panel.grid.minor.x = element_line(color = 'grey80', size = .1)
      # )
    
    #
    # save figures
    #
    
    f = file.path('output', 'figures', parameter_interp_rep)
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(plalt,
           filename = file.path(f, paste(tar_name(), '.pdf', sep = '')),
           width = 8, height = 6, dpi = 'print')
    
  }, 
  pattern = map(multiple_start_reps), 
  deployment = 'worker', 
  memory = 'transient',
  storage = 'worker'
)