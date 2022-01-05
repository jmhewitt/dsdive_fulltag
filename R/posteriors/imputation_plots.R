imputation_plots = function(mcmc_sample_dir, output_dir, template_bins, 
                            tag_info) {
  
  mcmc_sample_dir = 'output/mcmc/nim_fit'
  output_dir = 'output/new_plots/multinomial_prop_reduced_classes'
  targets::tar_load(template_bins)
  targets::tar_load(tag_info)
  targets::tar_load(movement_classes)
  
  library(ggplot2)
  library(ggthemes)
  library(dplyr)
  
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
                         ncol = length(nim_pkg$data$depths))
  for(f in stage_sample_files) {
    stage_samples = readRDS(f)
    stage_summary = stage_summary + apply(stage_samples, 2, function(samples) {
      tabulate(bin = samples, nbins = nim_pkg$consts$n_stages)
    })
  }
  rownames(stage_summary) =  rownames(movement_classes$stage_defs)
  
  # normalize posteriors
  stage_summary = stage_summary / colSums(stage_summary)[1]
  
  # enrich with information for plotting
  df = data.frame(
    t(stage_summary), 
    # time = as.POSIXct(nim_pkg$data$times, origin = '1970-01-01 00:00.00 UTC'), 
    depth_bin = nim_pkg$data$depths,
    depth = template_bins$center[nim_pkg$data$depths],
    sample_stages = rownames(movement_classes$stage_defs)[
      stage_samples[nrow(stage_samples),]
    ]
  ) %>% mutate(
    depth_lower = depth - template_bins$halfwidth[depth_bin],
    depth_upper = depth + template_bins$halfwidth[depth_bin]
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
    df$time = 1:nrow(df)
    df$stage_mode = rownames(stage_summary)[apply(stage_summary, 2, which.max)]
    
    # # plot posterior mode for stage speed
    # pl = ggplot(df[inds,] %>% 
    #               mutate(mode_dir = gsub(pattern = '_(ascent|descent|free)',  
    #                                      replacement = '', 
    #                                      x = stage_mode)), 
    #             aes(x = time, y = depth, col = mode_dir)) + 
    #   # imputed trajectory as path
    #   geom_line(col = 'grey80') +
    #   # imputed trajectory as discrete depths
    #   geom_point() + 
    #   # deep depth threshold
    #   geom_hline(yintercept = 800, lty = 3, alpha = .6) + 
    #   # formatting
    #   scale_color_brewer(type = 'qual', palette = 'Dark2') +
    #   # scale_color_manual(
    #   #   values = c(deep_descent = '#1f78b4', deep_ascent = '#33a02c',
    #   #              shallow_descent = '#a6cee3', shallow_ascent = '#b2df8a',
    #   #              deep_forage = '#d95f02', free_surface = '#e7298a')
    #   # ) + 
    #   depth_scale + 
    #   # scale_x_datetime(date_breaks = '12 hours', 
    #   #                  date_labels = c('%b %d', ' ')) + 
    #   theme_few() + 
    #   theme(axis.title.x = element_blank(), 
    #         legend.title = element_blank(),
    #         panel.grid.minor.y = element_line(color = 'grey90'),
    #         panel.grid.major.y = element_line(color = 'grey90'))
    
    # plot posterior mode for stage speed
    pl = ggplot(df[inds,] %>% 
                  mutate(mode_speed = gsub(pattern = '_(ascent|descent|free)',  
                                         replacement = '', 
                                         x = stage_mode),
                         mode_dir = gsub(pattern = '(fast|slow)_',  
                                         replacement = '', 
                                         x = stage_mode),
                         dir_angle = plyr::mapvalues(
                           x = mode_dir, 
                           from = c('ascent', 'descent', 'free'), 
                                    to = c(-45, -180+45, -90)),
                         speed_label = plyr::mapvalues(
                           x = mode_speed,
                           from = c('slow', 'fast'),
                           to = c('\u005E', '\u005E')# '\u2191')
                         ),
                         speed_face = plyr::mapvalues(
                           x = mode_speed,
                           from = c('slow', 'fast'),
                           to = c('plain', 'bold')
                         )
                        ), 
                aes(x = time, y = depth, ymin = depth_lower, ymax = depth_upper,
                    col = mode_speed, group = 1, angle = dir_angle, 
                    label = speed_label)) + 
      # imputed trajectory as path
      geom_line(alpha = .1, lwd = 1) +
      # imputed trajectory as discrete depths
      geom_point(pch = 3, alpha = .1) + 
      geom_linerange(alpha  = .1, lwd = 1) +
      geom_text() + 
      # deep depth threshold
      geom_hline(yintercept = 800, lty = 3, alpha = .6) + 
      # formatting
      scale_color_brewer(type = 'qual', palette = 'Dark2') +
      # scale_color_manual(
      #   values = c(deep_descent = '#1f78b4', deep_ascent = '#33a02c',
      #              shallow_descent = '#a6cee3', shallow_ascent = '#b2df8a',
      #              deep_forage = '#d95f02', free_surface = '#e7298a')
      # ) + 
      scale_color_manual(
        # values = c('fast' = '#ca0020', 'slow' = '#0571b0')
        # values = c('slow' = '#7a0177', 'fast' = '#fbb4b9')
        values = c(#'fast' = rgb(255/255, 69/255, 0/255), 
                   'fast' = rgb(220/255, 0/255, 5/255), 
                   # 'slow' = rgb(0,0,1),
                   'slow' = rgb(0/255,69/255,255/255))
      ) + 
      depth_scale + 
      # scale_x_datetime(date_breaks = '12 hours', 
      #                  date_labels = c('%b %d', ' ')) + 
      theme_few() + 
      theme(axis.title.x = element_blank(), 
            legend.title = element_blank())
    
    pl
    
    
    # save plot of dive record 
    f =  file.path(output_dir, 
                   paste('posterior_mode_speed_', subj_label, '.pdf', sep = '')
    )
    

    ggsave(pl, filename = f, dpi = 'print', height = 12, width = 12, 
           limitsize = FALSE)
    
    
    
    
    
    
    
    
    
    
    
    # plot posterior mode for stage direction
    pl = ggplot(df[inds,] %>% 
                  mutate(mode_dir = gsub(pattern = '(slow|medium|fast)_',  
                                           replacement = '', 
                                           x = stage_mode)), 
                aes(x = time, y = depth, col = mode_dir)) + 
      # imputed trajectory as path
      geom_line(col = 'grey80') +
      # imputed trajectory as discrete depths
      geom_point() + 
      # deep depth threshold
      geom_hline(yintercept = 800, lty = 3, alpha = .6) + 
      # formatting
      scale_color_brewer(type = 'qual', palette = 'Dark2') +
      # scale_color_manual(
      #   values = c(deep_descent = '#1f78b4', deep_ascent = '#33a02c',
      #              shallow_descent = '#a6cee3', shallow_ascent = '#b2df8a',
      #              deep_forage = '#d95f02', free_surface = '#e7298a')
      # ) + 
      depth_scale + 
      # scale_x_datetime(date_breaks = '12 hours', 
      #                  date_labels = c('%b %d', ' ')) + 
      theme_few() + 
      theme(axis.title.x = element_blank(), 
            legend.title = element_blank(),
            panel.grid.minor.y = element_line(color = 'grey90'),
            panel.grid.major.y = element_line(color = 'grey90'))
    
    # save plot of dive record 
    f =  file.path(output_dir, 
                   paste('posterior_mode_dir_', subj_label, '.pdf', sep = '')
    )
    
    ggsave(pl, filename = f, dpi = 'print', height = 12, width = 12*12, 
           limitsize = FALSE)
    
    
    # plot a sample stage imputation
    pl = ggplot(df[inds,] %>% 
                  mutate(sample_dir = gsub(pattern = '(slow|medium|fast)_',  
                                           replacement = '', 
                                           x = sample_stages)), 
                aes(x = time, y = depth, col = sample_dir)) + 
      # imputed trajectory as path
      geom_line(col = 'grey80') +
      # imputed trajectory as discrete depths
      geom_point() + 
      # deep depth threshold
      geom_hline(yintercept = 800, lty = 3, alpha = .6) + 
      # formatting
      scale_color_brewer(type = 'qual', palette = 'Dark2') +
      # scale_color_manual(
      #   values = c(deep_descent = '#1f78b4', deep_ascent = '#33a02c',
      #              shallow_descent = '#a6cee3', shallow_ascent = '#b2df8a',
      #              deep_forage = '#d95f02', free_surface = '#e7298a')
      # ) + 
      depth_scale + 
      # scale_x_datetime(date_breaks = '12 hours', 
      #                  date_labels = c('%b %d', ' ')) + 
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
    
    
  }
  
}
