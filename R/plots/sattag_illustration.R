sattag_illustration_script = tar_target(
  name = sattag_illustration,
  command = {
    
    #
    # illustration of raw sattag data
    #
    
    tag = 'ZcTag093'
    
    tar_load(template_bins)
    
    # clean up template bins for plotting
    for(i in 2:nrow(template_bins)) {
      # halfwidths must be non-decreasing
      template_bins$halfwidth[i] = 
        max(
          template_bins$halfwidth[i-1],
          template_bins$center[i] - 
            (template_bins$center[i-1] + template_bins$halfwidth[i-1])
        )
      # adjust bin center according to updated halfwidth
      template_bins$center[i] = 
        template_bins$center[i-1] + 
        template_bins$halfwidth[i-1] + 
        template_bins$halfwidth[i]
    }
  
    tag_ind = which(sapply(raw_sattags, function(x) x$tag) == tag)
    
    record = raw_sattags[[tag_ind]]
    
    # crudely compute covariates (i.e., without accounting for breaks, etc.)
    covariates = covariate_tx(
      covariates = rbind(
        daytime = daytime(date = record$times, lat = cape_hatteras_loc['lat'], 
                          lon = cape_hatteras_loc['lon']),
        moonlit = moonlit(date = record$times, lat = cape_hatteras_loc['lat'], 
                          lon = cape_hatteras_loc['lon']),
        time = record$times,
        depth = record$depths
      ), 
      control = covariate_tx_control
    )
    
    d = data.frame(
      # raw series data
      Date = record$times,
      bin = record$depth.bin,
      obs = 1:length(record$times),
      # add computed covariates
      t(covariates),
      # add raw celestial information (altitudes converted to degrees)
      moon = getMoonPosition(
        date = record$times,
        lat = cape_hatteras_loc['lat'],
        lon = cape_hatteras_loc['lon'],
        keep = 'altitude'
      )[,'altitude'] * 180 / pi,
      moon_illumination = getMoonIllumination(
        date = record$times,
        keep = 'fraction'
      )[,'fraction'],
      sun = getSunlightPosition(
        date = record$times,
        lat = cape_hatteras_loc['lat'],
        lon = cape_hatteras_loc['lon'],
        keep = 'altitude'
      )[,'altitude'] * 180 / pi
    ) %>% mutate(
      # enrich series data with bin centers and limits
      depth = template_bins$center[bin],
      lwr = depth - template_bins$halfwidth[bin],
      upr = depth + template_bins$halfwidth[bin]
    ) %>% filter(
      # focus the plotting window to specific time range
      record$exposure_time - duration(1.5, 'hour') <= Date,
      Date <= record$exposure_time + duration(1.5, 'hour')
    )
    
    # create a plottable celestial covariate
    d$celestial = (apply(d, 1, function(r) {
      unname(
        ifelse(
          as.numeric(r['daytime']) == TRUE, 
          'Daytime', 
          ifelse(as.numeric(r['night_dark']) == TRUE, 
                 'Dark night', 
                 'Moonlit night')
        )
      )
    }))
    
    model_times = c(d$Date, record$exposure_time)
    model_time_labels = parse(
      text = c(paste('t[', d$obs, ']', sep =''), 't~"*"')
    )
    o = order(model_times)
    
    pl = ggplot(d, aes(x = Date, y = depth, ymin = lwr, ymax = upr)) + 
      # horizontal surface line
      geom_hline(yintercept = template_bins$center[1], lty = 3) +
      # horizontal deep depth line
      geom_hline(yintercept = covariate_tx_control$deep_depth, 
                 lty = 2, alpha = .5) +
      # depth bin enumeration
      geom_hline(
        mapping = aes(yintercept = D),
        data = data.frame(
          D = c(
            (template_bins$center - template_bins$halfwidth)[1],
            template_bins$center + template_bins$halfwidth
          )
        ),
        alpha = .05
      ) +
      # time series plot of depth bins
      geom_line(col = 'grey60') +
      geom_pointrange() + 
      # cee start time
      geom_vline(xintercept = record$exposure_time, lty = 2, alpha = .5) + 
      # axis and plot formatting
      scale_x_datetime(
        'Time (UTC)', 
        date_breaks = 'hour',
        date_labels = '%H%M'#, 
        # # observation and CEE times as second axis
        # sec.axis = sec_axis(
        #   trans = identity,
        #   breaks = model_times[o], 
        #   labels = model_time_labels[o]
        #   
        # )
      ) + 
      scale_y_reverse(
        breaks = seq(0, max(template_bins$center + template_bins$halfwidth), 
                     by = 400), 
        limits = c(1400,0),
        # depth bins as second y-axis
        sec.axis = sec_axis(
          trans = identity,
          breaks = c(
            (template_bins$center - template_bins$halfwidth)[1],
            template_bins$center + template_bins$halfwidth
          ), 
          labels = parse(text = paste('italic(D)[', 0:16, ']', sep=''))
        )
      ) +
      ylab('Depth (m)') + 
      theme_few()
    
    # pl
    
    
    #
    # illustration of covariates
    #
    
    pl_vertical = ggplot(d, aes(x = Date, y = total_vertical_poly_1)) + 
      # time series plot of covariate
      geom_line(col = 'grey60') +
      geom_point() +
      # cee start time
      geom_vline(xintercept = record$exposure_time, lty = 2, alpha = .5) + 
      # axis and plot formatting
      scale_x_datetime(
        'Time (UTC)', 
        date_breaks = 'hour',
        date_labels = '%H%M'#, 
        # # observation and CEE times as second axis
        # sec.axis = sec_axis(
        #   trans = identity,
        #   breaks = model_times[o], 
        #   labels = model_time_labels[o]
        #   
        # )
      ) + 
      ylab('Recent vertical distance (m)') + 
      theme_few()
    
    # pl_vertical
    
    pl_prop_deep = ggplot(d, aes(x = Date, y = prop_deep_poly_1)) + 
      # time series plot of covariate
      geom_line(col = 'grey60') +
      geom_point() +
      # cee start time
      geom_vline(xintercept = record$exposure_time, lty = 2, alpha = .5) + 
      # axis and plot formatting
      scale_x_datetime(
        'Time (UTC)', 
        date_breaks = 'hour',
        date_labels = '%H%M'#, 
        # # observation and CEE times as second axis
        # sec.axis = sec_axis(
        #   trans = identity,
        #   breaks = model_times[o], 
        #   labels = model_time_labels[o]
        #   
        # )
      ) + 
      scale_y_continuous(
        'Proportion',
        limits = c(0,1), 
        breaks = round(
          x = seq(from = 0, to = 1, by = 3 * covariate_tx_control$obs_freq / 
                    covariate_tx_control$window_len), 
          digits = 2
        )
      ) + 
      theme_few()
    
    # pl_prop_deep
    
    pl_celestial = ggplot(
      data = d %>% 
        select(Date, sun, moon) %>% 
        pivot_longer(cols = sun:moon, names_to = 'series', 
                     values_to = 'altitude'),
      mapping = aes(x = Date, y = altitude, col = factor(series))
    ) + 
      # time series plot of covariate
      geom_line() +
      # cee start time
      geom_vline(xintercept = record$exposure_time, lty = 2, alpha = .5) + 
      # axis and plot formatting
      scale_x_datetime(
        'Time (UTC)', 
        date_breaks = 'hour',
        date_labels = '%H%M'#, 
        # # observation and CEE times as second axis
        # sec.axis = sec_axis(
        #   trans = identity,
        #   breaks = model_times[o], 
        #   labels = model_time_labels[o]
        #   
        # )
      ) + 
      ylab('Altitude (deg.)') +
      theme_few()
    
    # pl_celestial
    
    #
    # save output
    #
    
    title_size = 12
    
    pl_out = egg::ggarrange(
      pl + 
        theme(axis.title.x = element_blank(), 
              plot.title = element_text(hjust = .5, size = title_size)) + 
        ggtitle('Observed depth series') + 
        ylab('Depth (m)'), 
      pl_prop_deep + 
        theme(axis.title.x = element_blank(),
              plot.title = element_text(hjust = .5, size = title_size)) + 
        ggtitle(
          'Derived covariate: Proportion of last hour spent in deep depths'
        ),
      pl_vertical + 
        theme(plot.title = element_text(hjust = .5, size = title_size)) + 
        ggtitle(
          'Derived covariate: Total vertical distance traveled in last hour'
        ) + 
        ylab('Distance (m)'),
      ncol = 1, 
      heights = c(5,1,1)
    )
    
    # pl_out
    
    f = file.path('output','figures')
    dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
    ggsave(pl_out, filename = file.path(f, paste(tar_name(), '_', tag, '.pdf', 
                                             sep ='')),
           width = 8, height = 8, dpi = 'print')
  }
)