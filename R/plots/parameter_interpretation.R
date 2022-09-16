parameter_interpretation_plot_script = tar_target(
  name = parameter_interpretation_plot,
  command = {
    
    # parameter_interp_rep = paste(
    #   'parameter_interpretation_', multiple_start_reps, sep =''
    # )
    
    parameter_interp_rep = 'fixed_init_beta'
    
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
        t(rowMeans(samples$baseline_deep_pred_samples)),
        # add IDDI-like measure; contains end+start of deep dives + IDDI
        recovery = mean(apply(samples$baseline_deep_pred_samples, 2, 
                              function(r) {
          diff(r)
        }))
      )
    }))
    
    # TODO: temporarily identify missing simulations
    if(nrow(summaries) != 45) {
      message(paste('Missing simulation output for rep', multiple_start_reps))
    }
    
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
    summaries$second_deep = summaries$second_deep / 60
    summaries$recovery = summaries$recovery / 60
    
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
    
    # visualize parameter effects
    pl = ggplot(summaries, aes(x = time, y = first_deep, col = init_stage, 
                          lty = init_stage, pch = init_stage,
                          group = init_stage)) + 
      geom_point() + 
      geom_line() + 
      scale_linetype_discrete('Initial movement', labels = titlecase) + 
      scale_shape_discrete('Initial movement', labels = titlecase) + 
      scale_color_brewer('Initial movement', labels = titlecase, 
                         type = 'qual', palette = 'Dark2') + 
      xlab('Lighting conditions') + 
      ylab('Deep depth hitting time (min)') + 
      facet_wrap(~prototype, labeller = titlecase) + 
      theme_few()
    
    # alternate visualization for parameter effects
    plalt = ggplot(summaries, aes(x = prototype, y = first_deep, col = init_stage, 
                               lty = init_stage, pch = init_stage,
                               group = init_stage)) + 
      geom_point() + 
      geom_line() + 
      scale_linetype_discrete('Initial movement', labels = titlecase) + 
      scale_shape_discrete('Initial movement', labels = titlecase) + 
      scale_color_brewer('Initial movement', labels = titlecase, 
                         type = 'qual', palette = 'Dark2') + 
      xlab('Recent diving activity') + 
      # ylab('Deep depth hitting time (min)') + 
      ylab(expression(E~'['~H[1](1)~']')) + 
      facet_wrap(~time, labeller = titlecase) + 
      theme_few() + 
      theme(axis.title.y = element_text(angle = 0, vjust = .5))

    #
    # load and plot prototypes as well, since this helps interpret the panels
    #
    
    # load prototypes 
    ptypes = readRDS(
      file.path('output', 'parameter_interpretation', 
                'parameter_interpretation_patterns.rds')
    )
    
    # format prototypes for plotting
    df = do.call(rbind, apply(ptypes, 1, function(r) {
      data.frame(
        depth = r,
        t = ((seq_along(r) - 1) * covariate_tx_control$obs_freq)/60
      )
    })) %>% mutate(
      # add series labels
      ptype = factor(
        x = rev(rep(unlist(titlecase(levels(summaries$prototype))), 
                rep(ncol(ptypes), nrow(ptypes)))), 
        levels = unlist(titlecase(levels(summaries$prototype)))
      ),
      # map to depth bins
      bin = sapply(depth, function(d) {
        which.min(abs(d - template_bins$center))
      }),
      # add bin limits
      ctr = template_bins$center[bin],
      lwr = ctr - template_bins$halfwidth[bin],
      upr = ctr + template_bins$halfwidth[bin]
    )
    
    # visualize prototypes
    pl2 = ggplot(df, aes(x = t, y = ctr, ymin = lwr, ymax = upr)) + 
      # horizontal surface line
      geom_hline(yintercept = template_bins$center[1], lty = 3) +
      # horizontal deep depth line
      geom_hline(yintercept = covariate_tx_control$deep_depth, 
                 lty = 2, alpha = .5) +
      # depth trajectory
      geom_line(col = 'grey60') + 
      geom_pointrange() + 
      # panels
      facet_wrap(~ptype, nrow = 1) + 
      # formatting
      scale_y_reverse(
        'Depth (m)', 
        breaks = seq(from = 0, to = 1600, by = 400)
      ) + 
      scale_x_continuous(
        'Time (min)', 
        breaks = seq(from = 0, to = 60, by = 15)
      ) + 
      theme_few()
    
    
    #
    # load and plot summaries for alternate prototypes
    #
    
    alt_summaries = rbind(
      # summaries for alternates
      cbind(
        readRDS(
          file.path('output', 'parameter_interpretation', 
                    'parameter_interpretation_patterns_alt_summaries.rds')
        )[, c('prop_deep', 'total_vert')],
        chosen = FALSE
      ),
      # summaries for chosen
      t(apply(ptypes, 1, function(r) {
        c(prop_deep = mean(r >= 800), total_vert = sum(abs(diff(r))), 
          chosen = TRUE)
      }))
    )
    
    pl_summaries = ggplot(alt_summaries, aes(x = prop_deep, y = total_vert,
                                             col = factor(chosen))) + 
      geom_point() + 
      scale_color_manual(values = c('1' = '#d95f02', '0' = 'black')) + 
      xlab('Proportion deep') + 
      ylab('Total vertical distance') + 
      guides(color = 'none') + 
      theme_few()
      
    #
    # save figures
    #
    
    f = file.path('output', 'figures', parameter_interp_rep)
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    
    ggsave(plalt + theme(axis.text.x = element_text(size = 8)), 
           filename = file.path(f, paste(tar_name(), '.pdf', sep = '')),
           width = 8, height = 3, dpi = 'print')
    
    ggsave(pl2, filename = file.path(f, paste(tar_name(), '_prototypes.pdf', 
                                              sep = '')),
           width = 8, height = 4, dpi = 'print')
    
    ggsave(
      pl_summaries, 
      filename = file.path(f, paste(tar_name(), '_summaries.pdf', sep = '')),
      dpi = 'print', width = 6, height = 4
    )
    
    # output summary information for rep so we can compare findings
    summaries$rep = multiple_start_reps
    summaries
    
  }, 
  # pattern = map(multiple_start_reps), 
  deployment = 'worker', 
  memory = 'transient',
  storage = 'worker'
)