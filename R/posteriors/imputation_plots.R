imputation_plots = function(mcmc_sample_dir, output_dir, template_bins, 
                            tag_info) {
  
  # load source data
  nim_pkg = readRDS(dir(
    path = mcmc_sample_dir, pattern = 'nim_pkg.rds', full.names = TRUE
  ))
  
  # identify stage sample files
  stage_sample_files = dir(
    path = mcmc_sample_dir, pattern = 'stage.*.rds', full.names = TRUE
  )
  
  # skip plots if stages were not sampled
  if(length(stage_sample_files) == 0) {
    return(0)
  }
  
  # aggregate posterior stage samples
  stage_summary = matrix(as.integer(0), 
                         nrow = nim_pkg$consts$n_stages, 
                         ncol = nim_pkg$consts$n_timepoints)
  for(f in stage_sample_files) {
    stage_samples = readRDS(f)
    stage_summary = stage_summary + apply(stage_samples, 2, function(samples) {
      tabulate(bin = samples, nbins = nim_pkg$consts$n_stages)
    })
  }
  rownames(stage_summary) = names(nim_pkg$consts$movement_types)
  
  # normalize posteriors
  stage_summary = stage_summary / colSums(stage_summary)[1]
  
  # enrich with information for plotting
  df = data.frame(
    t(stage_summary), 
    time = as.POSIXct(nim_pkg$data$times, origin = '1970-01-01 00:00.00 UTC'), 
    depth_bin = nim_pkg$data$depths,
    depth = template_bins$center[nim_pkg$data$depths],
    sample_stages = names(nim_pkg$consts$movement_types)[
      stage_samples[nrow(stage_samples),]
    ],
    n_stage_support = colSums(nim_pkg$data$stage_supports)
  )
  
  # create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  depth_scale = scale_y_reverse('Depth (m)', minor_breaks = function(lim) {
    seq(from = round(min(lim)/50)*50, to = max(lim), by = 50)
  })
  
  # make plots for all tags
  for(subj_id in 1:nim_pkg$consts$n_subjects) {
    
    subj_label = nim_pkg$consts$subject_id_labels[subj_id]
    
    cee_start = tag_info %>% 
      dplyr::filter(deployid == subj_label) %>% 
      select(cee_start) %>% 
      unlist()
    
    # identify tag indices
    inds = apply(
      X = data.frame(nim_pkg$consts$segments) %>% 
        dplyr::filter(subject_id == subj_id), 
      MARGIN = 1, 
      FUN = function(r) {
        r['start_ind']:r['end_ind']
      }
    )
    
    if(is.list(inds)) {
      inds = do.call(c, inds)
    }
    
    # plot posterior probabilities for deep-forage state
    pl = ggplot(df[inds,],
                aes(x = time, y = depth, 
                    col = cut(deep_forage, 
                              breaks = seq(from = 0, to = 1, by = .2), 
                              include.lowest = TRUE))) + 
      # depth series
      geom_line(col = 'grey80') + 
      geom_point() +
      # deep depth indicator
      geom_hline(yintercept = 800, lty = 3) + 
      # cee time 
      geom_vline(xintercept = cee_start, lty = 3) + 
      # formatting
      depth_scale + 
      scale_x_datetime() + 
      scale_color_brewer('Deep-forage prob.', type = 'seq', palette = 'OrRd') + 
      theme_few() + 
      theme(axis.title.x = element_blank(),
            panel.grid.minor.y = element_line(color = 'grey90'),
            panel.grid.major.y = element_line(color = 'grey90'))
    
    f =  file.path(output_dir, 
                   paste('deep_forage_probs_', subj_label, '.pdf', sep = '')
    )
    
    ggsave(pl, filename = f, dpi = 'print', height = 12, width = 12*12, 
           limitsize = FALSE)
    
    # plot a sample stage imputation
    pl = ggplot(df[inds,], aes(x = time, y = depth, col = sample_stages)) + 
      # imputed trajectory as path
      geom_line(col = 'grey80') +
      # imputed trajectory as discrete depths
      geom_point() + 
      # deep depth threshold
      geom_hline(yintercept = 800, lty = 3, alpha = .6) + 
      # cee time 
      geom_vline(xintercept = cee_start, lty = 3) + 
      # formatting
      # scale_color_brewer(type = 'qual', palette = 'Dark2') + 
      scale_color_manual(
        values = c(deep_descent = '#1f78b4', deep_ascent = '#33a02c',
                   shallow_descent = '#a6cee3', shallow_ascent = '#b2df8a',
                   deep_forage = '#d95f02', free_surface = '#e7298a')
      ) + 
      depth_scale + 
      scale_x_datetime(date_breaks = '12 hours', 
                       date_labels = c('%b %d', ' ')) + 
      theme_few() + 
      theme(axis.title.x = element_blank(), 
            legend.title = element_blank(),
            panel.grid.minor.y = element_line(color = 'grey90'),
            panel.grid.major.y = element_line(color = 'grey90'))
    
    # save plot of dive record 
    f =  file.path(output_dir, 
                   paste('sample_imputation_', subj_label, '.pdf', sep = '')
    )
    
    ggsave(pl, filename = f, dpi = 'print', height = 12, width = 12*12, 
           limitsize = FALSE)
    
    # plot number of possible stages at each observation
    pl = ggplot(df[inds,], aes(x = time, y = depth, 
                               col = factor(n_stage_support))) + 
      # imputed trajectory as path
      geom_line(col = 'grey80') +
      # imputed trajectory as discrete depths
      geom_point() + 
      # deep depth threshold
      geom_hline(yintercept = 800, lty = 3, alpha = .6) + 
      # cee time 
      geom_vline(xintercept = cee_start, lty = 3) + 
      # formatting
      scale_color_viridis_d(direction = -1) + 
      depth_scale + 
      scale_x_datetime(date_breaks = '12 hours', 
                       date_labels = c('%b %d', ' ')) + 
      theme_few() + 
      theme(axis.title.x = element_blank(), 
            legend.title = element_blank(),
            panel.grid.minor.y = element_line(color = 'grey90'),
            panel.grid.major.y = element_line(color = 'grey90'))
    
    # save plot of dive record 
    f =  file.path(output_dir, 
                   paste('n_stages_possible_', subj_label, '.pdf', sep = '')
    )
    
    ggsave(pl, filename = f, dpi = 'print', height = 12, width = 12*12, 
           limitsize = FALSE)
    
  }
  
}
