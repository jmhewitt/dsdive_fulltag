parameter_interpretation_pattern_script = tar_target(
  name = parameter_interpretation_patterns, 
  command = {
    
    # export trajectory segments used to initialize posterior predictive 
    # simulations for parameter interpretation.
    # 
    # the script relies on clustering, in which the code does not reproducibly
    # generate cluster labels.  so, the output of the method cannot be 
    # automatically reproduced even though the key outputs essentially can be.
    #
    # as such, this target should be run manually.
    
    # tar_load(fit_marginalized_model)
    
    fit_marginalized_model = list(
      samples = file.path('output', 'mcmc', 'fit_marginalized_model'),
      package = file.path('output', 'mcmc', 'fit_marginalized_model', 
                          'nim_pkg.rds')
    )
    
    pkg = readRDS(fit_marginalized_model$package)
    
    
    #
    # study variety in diving trajectories 
    #
    
    # focus on times when animals were observed near surface
    surface_inds = which(pkg$data$depth_bins == 1)
    
    # extract covariates for these times
    surface_conditions = pkg$data$covariates[,surface_inds]
    
    # only retain times/data for which all covariate data is defined
    surface_inds = surface_inds[complete.cases(t(surface_conditions))]
    surface_conditions = pkg$data$covariates[,surface_inds]
    
    # extract the trajectories that preceded the near-surface observations
    ninds = covariate_tx_control$window_len / covariate_tx_control$obs_freq
    prior_trajectories = do.call(rbind, lapply(surface_inds, function(i) {
      pkg$data$depth_bins[seq(to = i - 1, length.out = ninds)]
    }))
    
    # basic clustering for the trajectory segments 
    clusters = kmeans(x = prior_trajectories, centers = 25, nstart = 10)
    
    # compute summary covariates for each cluster center
    covs = data.frame(
      t(apply(clusters$centers, 1, function(pattern) {
        depths = template_bins$center[round(pattern)]
        c('prop_deep' = mean(depths > covariate_tx_control$deep_depth),
          'total_vert' = sum(abs(diff(depths))))
      })), 
      # associate each cluster's number with the covariates
      pattern = 1:nrow(clusters$centers)
    )
    
    # # quickly explore the data
    # covs %>% arrange(prop_deep, total_vert, pattern)
    
    # select trajectories to serve as prototypes for interpreting parameters
    #  - the cluster id's are not replicable, but patterns essentially are
    patterns = c(
      4,  # has been doing shallow diving for a while
      20, # finishing a deep dive, and starting recovery
      2   # immediately coming out of a deep dive
    )
    
    # breakpoint at which we want to manually check the quality of the patterns
    # due to reproducibility issues in the labels
    browser()
    
    # compare summaries of selected trajectories against all cluster centers
    plot(total_vert ~ prop_deep, covs)
    points(covs[patterns,'prop_deep'], covs[patterns,'total_vert'], col = 2,
           pch = 16)
    
    # plottable form for clustered trajectories
    df = data.frame(round(clusters$centers)) %>% 
      pivot_longer(
        cols = everything(), 
        names_to = 'window_ind', 
        values_to = 'depth_bin'
      ) %>% mutate(
        window_ind = as.numeric(
          gsub(pattern = 'X', replacement = '', x = window_ind)
        ),
        depth = template_bins$center[depth_bin],
        id = rep(1:nrow(clusters$centers), 
                 rep(ncol(clusters$centers), nrow(clusters$centers)))
      )
    
    # add depth bin limits
    df = df %>% mutate(
      lwr = depth - template_bins$halfwidth[depth_bin],
      upr = depth + template_bins$halfwidth[depth_bin],
      t = (window_ind - 1) * 5
    )
      
    # plot time series for each cluster center
    pl = ggplot(df, aes(x = t, y = depth, ymin = lwr, ymax = upr)) + 
      # horizontal surface line
      geom_hline(yintercept = template_bins$center[1], lty = 3) +
      # horizontal deep depth line
      geom_hline(yintercept = covariate_tx_control$deep_depth, 
                 lty = 2, alpha = .5) +
      # depth trajectory
      geom_line(col = 'grey60') + 
      geom_pointrange(size = .3) + 
      # panels
      facet_wrap(~factor(id), ncol = 5) + 
      # formatting
      scale_y_reverse(
        'Depth (m)', 
        breaks = seq(from = 0, to = 1600, by = 400)
      ) + 
      scale_x_continuous(
        'Time (min)', 
        breaks = seq(from = 0, to = 60, by = 15), 
        limits = c(0, 60)
      ) + 
      theme_few()
    
    pl
    
    # format prototypes for saving
    ptypes = round(clusters$centers)[patterns,]
    ptypes = matrix(template_bins$center[ptypes], nrow = nrow(ptypes))
    
    # recall that all of these end at the surface
    ptypes = cbind(ptypes, template_bins$center[1])
    
    # label entries
    colnames(ptypes) = paste('t', 1:ncol(ptypes), sep ='')
    rownames(ptypes) = paste('prototype', 1:nrow(ptypes), sep ='')
    
    # save prototype trajectories
    f = file.path('output', 'parameter_interpretation')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    saveRDS(ptypes, file = file.path(f, paste(tar_name(), '.rds', sep = '')))
    
    # save summaries for additional cluster centers
    saveRDS(covs, file = file.path(f, paste(tar_name(), '_alt_summaries.rds', 
                                            sep = '')))
    
    # save plots of additional cluster centers
    f = file.path('output', 'figures')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    ggsave(pl, filename = file.path(f, paste(tar_name(), '_alt.pdf', sep = '')), 
           width = 8, height = 10, dpi = 'print')
    
    
    
    f
  }, 
  # make sure this target does not run automatically; comment line out to run 
  # manually
  cue = tar_cue(mode = 'never')
)
