plot_targets = list(

  tar_target(
    name = tag_info_table,
    command = {
      do.call(rbind, lapply(raw_sattags, function(tag) { 
        data.frame(
          tag = tag$tag,
          date_start = date(min(tag$times)),
          date_end = date(max(tag$times)),
          cee_start = tag$exposure_time
        )
      }))
    }
  ),
  
  tar_target(
    name = zc093_figure, 
    command = {
      
      browser()
      # get index for zc093
      zc93ind = which(tag_info$deployid == 'ZcTag093')
      
      # window around cee, for plotting
      cee_start = tag_info$cee_start[zc93ind]
      cee_window = cee_start + duration(num = c(-6,1), units = 'hours')
      
      # depth data
      dat = raw_sattags[[zc93ind]]
      
      # depth profile around exposure time
      df = data.frame(
        time = dat$times,
        depth = template_bins$center[dat$depth.bin],
        depth_lwr = template_bins$center[dat$depth.bin] - 
          template_bins$halfwidth[dat$depth.bin],
        depth_upr = template_bins$center[dat$depth.bin] + 
          template_bins$halfwidth[dat$depth.bin]
      ) %>% 
        filter(cee_window[1] <= time, time <= cee_window[2])
      
      pl = ggplot(df, aes(x = time, ymin = depth_lwr, ymax = depth_upr, 
                          y = depth)) + 
        geom_hline(yintercept = 0, col = 'grey90') + 
        geom_line(col = 'grey80') + 
        geom_pointrange(size = .2) + 
        scale_y_reverse() + 
        scale_x_datetime(breaks = 'hour', date_labels = '%H:%M') + 
        geom_vline(xintercept = cee_start, lty = 3) + 
        xlab('Time (UTC)') + 
        ylab('Depth (m)') + 
        theme_few() + 
        theme()
      
      # create output directory and file name
      d = file.path('output', 'figures')
      dir.create(path = d, showWarnings = FALSE, recursive = TRUE)
      f = file.path('output', 'figures', 'zc93.png')
      
      # save file
      ggsave(pl, filename = f, dpi = 'print', width = 8, height = 3)
      
      pl
    }
  )
)
