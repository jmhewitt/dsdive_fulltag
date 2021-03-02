impute_observations = function(tag, endpoints, timestep, imputation_factor,
                               template_bins, stages, surface_run_length,
                               imputed_dive_label_plot_dir) {
  # Parameters:
  #  tag - data about satellite tags
  #  endpoints - information about the likely start/end times of dives
  #  timestep - number of seconds between observations
  #  imputation_factor - (integer) factor used to upscale sampling rate
  #  template_bins - reference depth bins
  #  surface_run_length - number of observations needed to label a run of 
  #   surface depth bin observations as a "free_surface" period
  
  if(imputation_factor <= 1) {
    stop(paste('Imputation factor must be larger than 1 to ensure stationary',
               'movement between observations.', sep = ' '))
  }
  
  # remove outermost list wrappings
  tag = tag[[1]]
  endpoints = endpoints[[1]]
  
  # initialize imputation with refined times and observed bins
  imputed = data.frame(
    time = seq(from = min(tag$times), to = max(tag$times), 
               by = timestep / imputation_factor)
  ) %>% dplyr::left_join(
    data.frame(time = tag$times, bin = tag$depth.bin, depth = tag$depths), 
    by = 'time'
  ) %>% dplyr::mutate(
    ind = 1:n()
  )
  
  # randomly assign dive starts/ends at discretized times (without replacment)
  endpoint_inds = numeric(nrow(endpoints$endpoint.inds))
  for(i in 1:length(endpoint_inds)) {
    r = endpoints$endpoint.inds[i,]
    x = imputed %>% 
      dplyr::filter(
        r['t_lwr'] <= time, time <= r['t_upr'], 
        any( is.na(bin), # allow transitions between observations
             bin == 1 ), # surface obs. can be transitions
        !(ind %in% endpoint_inds) # endpoints must be unique!
      ) %>% 
      dplyr::sample_n(1) %>%
      dplyr::select(ind) %>% 
      unlist()
    endpoint_inds[i] = x
  }
  
  # add the imputed surface observations and depths
  imputed$bin[endpoint_inds] = 1
  imputed$depth[endpoint_inds] = template_bins$center[1]
  
  # forward difference: time elapsed from one observation to the next
  time_increments = difftime(
    time2 = tag$times[1:(length(tag$times)-1)], 
    time1 = tag$times[2:length(tag$times)], 
    units = 'secs'
  )
  
  # linearly impute depth sequence
  df.observed = imputed %>% dplyr::filter(is.finite(depth))
  df.impute = imputed %>% dplyr::filter(is.na(depth))
  df.impute$depth = approx(
    x = df.observed$time, y = df.observed$depth, xout = df.impute$time
  )$y
  
  # map imputed depths to standardized bins
  df.impute$bin = sapply(df.impute$depth, function(depth) {
    which.min(abs(depth - template_bins$center))
  })
  df.impute$depth = template_bins$center[df.impute$bin]
  imputed[df.impute$ind, c('bin', 'depth')] = df.impute[, c('bin', 'depth')]
  imputed$observed = imputed$time %in% tag$times
  
  # identify (inclusive) time ranges for data gaps
  gap_after = which(time_increments > timestep)
  gap_ranges = data.frame(
    start = tag$times[gap_after] + timestep,
    end = tag$times[gap_after + 1] - timestep
  )
  
  # remove imputed depths during data gaps, and label gaps in data
  imputed$data_gap = FALSE
  if(nrow(gap_ranges) > 0) {
    for(i in 1:nrow(gap_ranges)) {
      gap_times = (gap_ranges[i, 'start'] <= imputed$time) & 
        (imputed$time <= gap_ranges[i, 'end'])
      imputed$bin[gap_times] = NA
      imputed$depth[gap_times] = NA
      imputed$data_gap[gap_times] = TRUE
    }
  }
  
  #
  # compute covariates for day/night and lunar illumination
  #

  # enrich data with celestial covariates
  imputed$daytime = daytime(
    date = imputed$time, 
    lat = tag$proximal_loc['lat'], 
    lon = tag$proximal_loc['lon']
  )
  imputed$moonlit = moonlit(
    date = imputed$time, 
    lat = tag$proximal_loc['lat'], 
    lon = tag$proximal_loc['lon']
  )

  
  #
  # supports for possible dive stage distributions
  #
  
  stage_support = matrix(0, nrow = length(stages), ncol = nrow(imputed))
  rownames(stage_support) = names(stages)
  
  for(i in 1:nrow(endpoints$dive.ranges)) {
    # indices for dive
    start_ind = endpoint_inds[endpoints$dive.ranges$T0_endpoint[i]]
    end_ind = endpoint_inds[endpoints$dive.ranges$T3_endpoint[i]]
    dive_inds = start_ind:end_ind
    # assign stages according to dive type
    if(endpoints$dive.ranges$type[i] == 'Shallow') {
      # identify observations before midpoint
      descent_inds = dive_inds < median(dive_inds)
      # disallow ascent movement before midpoint, also allow the possibility of 
      # coming out of a deep dive
      stage_support[
        c('deep_ascent', 'shallow_descent', 'free_surface'),
        dive_inds[descent_inds]
      ] = 1
      stage_support[
        c('deep_descent', 'shallow_ascent', 'free_surface'),
        dive_inds[!descent_inds]
      ] = 1
    } else {
      # extract depths
      dive_depths = imputed$depth[dive_inds]
      # find max depth, and indices surrounding max depth
      max.depth = max(dive_depths, na.rm = TRUE)
      foraging.range = range(which(dive_depths == max.depth))
      # assume foraging around the deepest portion of the dive
      forage_inds = dive_inds[
        seq(from = foraging.range[1], to = foraging.range[2])
      ]
      stage_support['deep_forage', forage_inds] = 1
      # disallow deep_ascent before foraging
      pre_forage_inds = dive_inds[dive_inds < min(forage_inds)]
      stage_support[
        c('deep_descent', 'deep_forage', 'shallow_descent', 'shallow_ascent', 
          'free_surface'), 
        pre_forage_inds
      ] = 1
      # disallow deep_descent after foraging
      post_forage_inds = dive_inds[dive_inds > max(forage_inds)]
      stage_support[
        c('deep_forage', 'deep_ascent', 'shallow_descent', 'shallow_ascent', 
          'free_surface'),
        post_forage_inds
      ] = 1
    }
  }
  
  # assign (remaining) obs not associated with dives to unrestricted movement
  stage_support[, colSums(stage_support) == 0] = 1
  
  # remove last observations before gap, and first observations after gap.
  # this step should be done after initially assigning stage supports so that
  # no stage support information is lost.  we only want to exclude data from 
  # estimation.
  if(nrow(gap_ranges) > 0) {
    for(i in 1:nrow(gap_ranges)) {
      # remove last observations before data gap
      remove_from = (endpoint_inds + 1)[max(
        which(imputed$time[endpoint_inds + 1] <= gap_ranges[i, 'start'])
      )]
      if(is.finite(remove_from)) {
        remove_inds = seq(
          from = remove_from,
          to = which(imputed$time == gap_ranges[i, 'start'])
        )
        imputed$bin[remove_inds] = NA
        imputed$depth[remove_inds] = NA
        imputed$data_gap[remove_inds] = TRUE
      }
      # remove first observations after data gap
      remove_to = (endpoint_inds - 1)[min(
        which(gap_ranges[i, 'end'] <= imputed$time[endpoint_inds - 1])
      )]
      if(is.finite(remove_to)) {
        remove_inds = seq(
          from = which(imputed$time == gap_ranges[i, 'end']),
          to = remove_to
        )
        imputed$bin[remove_inds] = NA
        imputed$depth[remove_inds] = NA
        imputed$data_gap[remove_inds] = TRUE
      }
    }
  }
  
  # 
  # label dive stages
  #
  
  imputed$stage = NA
  
  
  # label deep/shallow descent/forage/ascent dive stages
  for(i in 1:nrow(endpoints$dive.ranges)) {
    # indices for dive
    start_ind = endpoint_inds[endpoints$dive.ranges$T0_endpoint[i]]
    end_ind = endpoint_inds[endpoints$dive.ranges$T3_endpoint[i]]
    dive_inds = start_ind:end_ind
    # assign stages according to dive type
    if(endpoints$dive.ranges$type[i] == 'Shallow') {
      # for shallow dives, assume descending before dive midpoint
      descent_inds = dive_inds < median(dive_inds)
      # assign stages
      imputed$stage[
        dive_inds[descent_inds]
      ] = stages['shallow_descent']
      imputed$stage[
        dive_inds[!descent_inds]
      ] = stages['shallow_ascent']
    } else {
      # extract depths
      dive_depths = imputed$depth[dive_inds]
      # find max depth, and stage threshold via 85% rule
      max.depth = max(dive_depths)
      stage.thresh = .85 * max.depth
      # compute observed stage vector
      bottom.range = range(which(dive_depths >= stage.thresh))
      descent_inds = dive_inds[dive_inds < dive_inds[bottom.range[1]]]
      ascent_inds = dive_inds[dive_inds > dive_inds[bottom.range[2]]]
      forage_inds = dive_inds[!(dive_inds %in% c(descent_inds, ascent_inds))]
      # assign stages
      imputed$stage[descent_inds] = stages['deep_descent']
      imputed$stage[ascent_inds] = stages['deep_ascent']
      imputed$stage[forage_inds] = stages['deep_forage']
    }
  }
  
  # assign (remaining) obs not associated with dives to unrestricted movement
  imputed$stage[
    (imputed$data_gap == FALSE) & is.na(imputed$stage)
  ] = stages['free_surface']
  
  # munge time class covariate information for plotting
  imputed$time_class = factor(unlist(apply(imputed, 1, function(r) {
    if(r['daytime'] == TRUE) {
      'Daylight'
    } else {
      if(r['moonlit'] == TRUE) {
        'Moonlit'
      } else {
        'Night'
      }
    }
  })))
  
  # plot tag with imputed labels
  pl = ggplot(imputed %>% 
                mutate(stage = factor(stage, labels = names(stages))), 
              aes(x = time, y = depth, col = stage)) + 
    # underlay time class bands
    geom_rect(mapping = aes(xmin = time - 30, xmax = time + 30, 
                            ymin = 0, 
                            ymax = sum(template_bins[nrow(template_bins),]),
                            fill = time_class),
              alpha = .09, inherit.aes = FALSE) + 
    # cee time
    geom_vline(xintercept = tag$exposure_time, lty = 3, alpha = .6) + 
    # imputed trajectory as path
    geom_line(col = 'grey80') +
    # highlight observations
    # geom_point(data = imputed %>% dplyr::filter(observed == TRUE), 
    #            mapping = aes(x = time, y = depth), inherit.aes = FALSE,
    #            size = 1) + 
    # imputed trajectory as discrete depths
    geom_point() + 
    # deep depth threshold
    geom_hline(yintercept = 800, lty = 3, alpha = .6) + 
    # formatting 
    # scale_color_brewer(type = 'qual', palette = 'Dark2') + 
    scale_color_manual(
      values = c(deep_descent = '#1f78b4', deep_ascent = '#33a02c',
                 shallow_descent = '#a6cee3', shallow_ascent = '#b2df8a',
                 deep_forage = '#d95f02', free_surface = '#e7298a')
    ) + 
    scale_fill_manual(
      values = c(Daylight = '#ffc26c', Moonlit = '#75baff', Night = '#000000')
    ) + 
    scale_y_reverse('Depth (m)') + 
    scale_x_datetime(date_breaks = '12 hours', 
                     date_labels = c('%b %d', ' ')) + 
    theme_few() + 
    theme(axis.title.x = element_blank(), 
          legend.title = element_blank())
  
  # save plot of dive record 
  ggsave(pl, filename = file.path(imputed_dive_label_plot_dir, 
                                  paste(tag$tag, '.pdf', sep = '')),
         dpi = 'print', height = 12, width = 12*12, limitsize = FALSE)
  
  # package results
  list(
    tag = tag$tag,
    depth.bin = imputed$bin,
    depths = imputed$depth,
    times = imputed$time,
    daytime = imputed$daytime,
    moonlit = imputed$moonlit,
    stages = imputed$stage,
    stage_support = stage_support,
    data_gaps = imputed$data_gap,
    exposure_time = tag$exposure_time,
    baseline_end = tag$baseline_end
  )
}