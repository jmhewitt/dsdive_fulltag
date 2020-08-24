validation_statistics = function(nim_pkg, samples, tag, deep_dive_depth, 
                                 template_bins, sattag_timestep, 
                                 validation_output_dir, tag_sex) {

  if(file.exists('nim_pkg_0.5.rds')) {
    nim_pkg = readRDS('nim_pkg_0.5.rds')
  }
  
  # backwards compatability with existing code
  depth.bins = template_bins
  
  # create output directory
  o = file.path(validation_output_dir, tag)
  dir.create(o, showWarnings = FALSE, recursive = TRUE)
  
  # get tag id from target name
  tag_id = as.numeric(str_extract(tag, '[0-9]+'))
  
  # extract validation dives
  dives.obs = apply(
    X = data.frame(nim_pkg$consts$dive_relations_validation) %>% 
      dplyr::filter(tag == tag_id), 
    MARGIN = 1, 
    FUN = function(r) {
      inds = r['depth_first']:r['depth_last']
      list(
        depth.bins = template_bins,
        dive = list(depths = nim_pkg$data$depths[inds],
                    times = nim_pkg$data$times[inds] - 
                      nim_pkg$data$times[inds][1])
      )
    }
  )
  
  # format validation dives
  dives.obs.list = lapply(dives.obs, function(d) d$dive)
  
  # use 85% rule to approximate stage transition times (in minutes)
  times.stages.est = times.stages(dives.obs = dives.obs)
  
  # extract information about validation dives
  n.dives = length(dives.obs.list)
  validation.obs = do.call(rbind, lapply(1:n.dives, function(ind) {
    # extract dive
    d = dives.obs.list[[ind]]
    # extract estimates of stage durations
    stages.dur = times.stages.est[ind,]*60
    stages.dur[3] = d$times[length(d$times)] - sum(stages.dur)
    # loop over all observation times
    do.call(rbind, lapply(sattag_timestep, function(tstep) {
      # observed duration of dive
      duration.obs = d$times[length(d$times)]
      # extract summary features of dives
      data.frame(
        dive.id = ind,
        duration.obs = d$times[length(d$times)],
        max.depth.obs = template_bins$center[max(d$depths)],
        n.obs = length(d$times),
        n.tx.obs = sum(diff(d$depths) != 0),
        tstep = sattag_timestep,
        dur.s1 = as.numeric(stages.dur[1]),
        dur.s2 = as.numeric(stages.dur[2]),
        dur.s3 = as.numeric(stages.dur[3])
      )
    }))
  }))
  
  # compute depth bin distributions over time
  validation.depthseq = do.call(rbind, lapply(1:n.dives, function(ind) {
    # extract dive
    d = dives.obs.list[[ind]]
    # return tidy data frame with depth bins over time
    data.frame(dive.id = ind, depth = depth.bins$center[d$depths], 
               time = d$times)
  }))
  
  # extract information about posterior predictive samples
  n.samples = length(samples)
  postpred.samples.obs = do.call(rbind, lapply(1:n.samples, function(ind) {
    # extract dive
    d = samples[[ind]]
    # loop over all observation times
    do.call(rbind, lapply(sattag_timestep, function(tstep) {
      # exact duration of dive
      duration = d$dive$times[length(d$dive$times)] - d$dive$times[1]
      # extract summary features of dives
      data.frame(
        sample = ind,
        duration = duration,
        duration.obs = d$dive.obs$times[length(d$dive.obs$times)] - 
          d$dive.obs$times[1],
        depth.start.obs = d$dive.obs$depths[1],
        max.depth = depth.bins$center[max(d$dive$depths)],
        max.depth.obs = depth.bins$center[max(d$dive.obs$depths)],
        n.obs = length(d$dive.obs$times),
        n.tx = sum(diff(d$dive$depths) != 0),
        n.tx.obs = sum(diff(d$dive.obs$depths) != 0),
        tstep = tstep,
        dur.s1 = d$stages.dur[1],
        dur.s2 = d$stages.dur[2],
        dur.s3 = d$stages.dur[3]
      )
    }))
  }))
  
  # compute depth bin distributions over time
  postpred.depthseq = do.call(rbind, lapply(1:n.samples, function(ind) {
    # extract dive
    d = samples[[ind]]
    # return tidy data frame with depth bins over time
    data.frame(sample = ind,
               depth = depth.bins$center[d$dive.obs$depths],
               duration.obs = d$dive.obs$times[length(d$dive.obs$times)],
               depth.start.obs = d$dive.obs$depths[1],
               time = d$dive.obs$times,
               max.depth.obs = depth.bins$center[max(d$dive.obs$depths)])
  }))

  
  #
  # max observed depth distributions
  #
  
  df = rbind(
    postpred.samples.obs %>%
      dplyr::filter(max.depth.obs >= deep_dive_depth,
                    depth.start.obs == 1) %>%
      mutate(series = 'Post. Predictive',
             total = length(unique(sample)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      group_by(max.depth.obs, series) %>%
      summarise(prob = n() / total[1],
                n = n(),
                eps = eps[1]) %>%
      group_by(series) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1)),
    validation.obs %>%
      mutate(series = 'Empirical Validation',
             total = length(unique(dive.id)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      group_by(max.depth.obs, series) %>%
      summarise(prob = n() / total[1],
                n = n(),
                eps = eps[1]) %>%
      group_by(series) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1))
  )
  
  pl = ggplot(df, aes(x = max.depth.obs, y = cdf, ymin = cdf.lwr, 
                      ymax = cdf.upr, fill = series, col = series)) +
    geom_ribbon(alpha = .05, col = NA) +
    geom_point() +
    geom_line(lty = 3, alpha = .6) +
    scale_color_brewer('Distribution', type = 'qual', palette = 'Dark2') +
    scale_fill_brewer('Distribution', type = 'qual', palette = 'Dark2') +
    xlab('Max. observed depth (m)') +
    ylab('CDF') +
    theme_few() +
    theme(panel.border = element_blank())
  
  ggsave(pl, filename = file.path(o, 'max_observed_depth_cdf.png'),
         dpi = 'print')
  
  sink(file.path(o, paste('max_observed_depth_imse.txt')))
  r.imse = imse.gof(df, 'max.depth.obs')
  sink()
  
  sink(file.path(o, paste('max_observed_depth_cvm.txt')))
  r.imse = cvm.gof(df, 'max.depth.obs')
  sink()
  
  sink(file.path(o, paste('max_observed_depth_chisq.txt')))
  r = chisq.gof(df, 'max.depth.obs', collapse = TRUE)
  sink()
  
  sink(file.path(o, paste('max_observed_depth_chisq_collapsed.txt')))
  r = chisq.gof(df, 'max.depth.obs', collapse = TRUE)
  sink()
  
  sink(file.path(o, paste('max_observed_depth_ks.txt')))
  r = ks.gof(df, 'max.depth.obs')
  sink()
  
  save(df, r, r.imse, file = file.path(o, paste('max_observed_depth.RData')))
  
  
  #
  # observed duration distributions
  #
  
  df = rbind(
    postpred.samples.obs %>%
      dplyr::filter(max.depth.obs >= deep_dive_depth,
                    depth.start.obs == 1) %>%
      mutate(series = 'Post. Predictive',
             total = length(unique(sample)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      group_by(duration.obs, series) %>%
      summarise(prob = n() / total[1],
                eps = eps[1]) %>%
      group_by(series) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1)),
    validation.obs %>%
      mutate(series = 'Empirical Validation',
             total = length(unique(dive.id)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      group_by(duration.obs, series) %>%
      summarise(prob = n() / total[1],
                n = n(),
                eps = eps[1]) %>%
      group_by(series) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1))
  )
  
  pl = ggplot(df, aes(x = duration.obs/60, y = cdf, col = series, fill = series,
                      ymin = cdf.lwr, ymax = cdf.upr)) +
    geom_ribbon(alpha = .05, col = NA) +
    geom_point() +
    geom_line(lty = 3, alpha = .6) +
    scale_color_brewer('Distribution', type = 'qual', palette = 'Dark2') +
    scale_fill_brewer('Distribution', type = 'qual', palette = 'Dark2') +
    xlab('Observed dive duration (min)') +
    ylab('CDF') +
    theme_few() +
    theme(panel.border = element_blank())
  
  ggsave(pl, filename = file.path(o, 'observed_duration_cdf.png'),
         dpi = 'print')
  
  sink(file.path(o, paste('observed_duration_imse.txt')))
  r.imse = imse.gof(df, 'duration.obs')
  sink()
  
  sink(file.path(o, paste('observed_duration_chisq.txt')))
  r = chisq.gof(df, 'duration.obs', collapse = FALSE)
  sink()
  
  sink(file.path(o, paste('observed_duration_chisq_collapsed.txt')))
  r = chisq.gof(df, 'duration.obs', collapse = TRUE)
  sink()
  
  sink(file.path(o, paste('observed_duration_ks.txt')))
  r = ks.gof(df, 'duration.obs')
  sink()
  
  save(df, r, r.imse, file = file.path(o, paste('observed_duration.RData')))
  
  
  #
  #  distribution of depths at observation times
  #
  
  df = rbind(
    postpred.depthseq %>%
      dplyr::filter(max.depth.obs >= deep_dive_depth,
                    depth.start.obs == 1) %>%
      group_by(time) %>%
      mutate(total = length(unique(sample)),
             eps = sqrt(log(2/.05)/(2*total)),
             Distribution = 'Post. Predictive') %>%
      group_by(time, depth, Distribution) %>%
      summarise(prob = n() / total[1],
                n = n(),
                eps = eps[1]) %>%
      group_by(Distribution, time) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1)),
    validation.depthseq %>%
      group_by(time) %>%
      mutate(total = length(unique(dive.id)),
             eps = sqrt(log(2/.05)/(2*total)),
             Distribution = 'Empirical Validation') %>%
      group_by(time, depth, Distribution) %>%
      summarise(prob = n() / total[1],
                n = n(),
                eps = eps[1]) %>%
      group_by(Distribution, time) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1))
  )
  
  pl = ggplot(df, aes(x = (depth), y = cdf, ymin = cdf.lwr, ymax = cdf.upr,
                      fill = Distribution, col = Distribution,
                      group = Distribution)) +
    geom_ribbon(alpha = .05, col = NA) +
    geom_point() +
    geom_line(lty = 3, alpha = .6) +
    scale_color_brewer(type = 'qual', palette = 'Dark2') +
    scale_fill_brewer(type = 'qual', palette = 'Dark2') +
    xlab('Depth (m)') +
    ylab('CDF') +
    facet_wrap(~factor(time/60)) +
    # labels = paste(unique(time)/60, 'min', sep = ' '))) +
    theme_few()
  
  sink(file.path(o, paste('depths_by_time_ks.txt')))
  r = lapply(sort(unique(df$time))[-1], function(s) {
    tryCatch({
      res = ks.gof(df %>% mutate(series = Distribution) %>%
                     dplyr::filter(time==s), 'depth')
      cat(paste("\n(t=", s/60, " min. results)\n\n\n", sep=''))
      res
    }, error = function(e){ cat("ERROR :",conditionMessage(e), "\n") })
  })
  sink()
  
  sink(file.path(o, paste('depths_by_time_cvm.txt')))
  r = lapply(sort(unique(df$time))[-1], function(s) {
    tryCatch({
      res = cvm.gof(df %>% mutate(series = Distribution) %>%
                     dplyr::filter(time==s), 'depth')
      cat(paste("\n(t=", s/60, " min. results)\n\n\n", sep=''))
      res
    }, error = function(e){ cat("ERROR :",conditionMessage(e), "\n") })
  })
  sink()
  
  pvals = sapply(r[1:12], function(res) {
    ifelse(res$p.value < .01, '(p < .01)', 
           paste('(p = ', round(res$p.value, 2), ')', sep =''))
  })
  
  pvals.numeric = sapply(r[1:12], function(res) {
    res$p.value
  })
  
  time_labeller = function(labels, multi_line = TRUE)
  {
    labels <- lapply(labels, function(lab) {
      # paste('t =', lab, 'min.\n', pvals, sep = ' ')
      paste('t =', lab, 'min.', sep = ' ')
    })
    if (multi_line) {
      labels
    }
    else {
      collapse_labels_lines(labels)
    }
  }
  
  
  pl.pub = ggplot(df %>% dplyr::filter(time/60 >= 5, time/60 <= 60),
                  aes(x = (depth), y = cdf, ymin = cdf.lwr, ymax = cdf.upr,
                      fill = Distribution, col = Distribution,
                      group = Distribution)) +
    # geom_ribbon(alpha = .1, col = NA) +
    geom_point(size = 1) +
    geom_line(lty = 1, alpha = .2, lwd = .5) +
    scale_color_brewer(type = 'qual', palette = 'Dark2') +
    scale_fill_brewer(type = 'qual', palette = 'Dark2') +
    xlab('Depth (m)') +
    ylab(expression(P('\u2113'(t) <= x))) +
    facet_wrap(~factor(time/60), labeller = time_labeller) +
    theme_few()
  
  ggsave(pl.pub, filename = file.path(o, 'depths_by_time_cdf_pub.png'),
         dpi = 'print', width = 8, height = 6)
  
  ggsave(pl.pub, 
         filename = file.path(o, paste('depths_by_time_cdf_pub_', 
                                       tag_sex$deployid[tag_id], '.png', 
                                       sep = '')),
         dpi = 'print', width = 8, height = 6)
  
  ggsave(pl, filename = file.path(o, 'depths_by_time_cdf.png'),
         dpi = 'print', width = 14, height = 14)
  
  sink(file.path(o, paste('depths_by_time_imse.txt')))
  r.imse = lapply(sort(unique(df$time))[-1], function(s) {
    res = imse.gof(df %>% mutate(series = Distribution) %>%
                     dplyr::filter(time==s), 'depth')
    cat(paste("\n(t=", s/60, " min. results)\n\n\n", sep=''))
    c(s = s, imse = res)
  })
  sink()
  
  sink(file.path(o, paste('depths_by_time_chisq.txt')))
  r = lapply(sort(unique(df$time))[-1], function(s) {
    tryCatch({
      res = chisq.gof(df %>% mutate(series = Distribution) %>%
                        dplyr::filter(time==s), 'depth', collapse = FALSE)
      cat(paste("\n(t=", s/60, " min. results)\n\n\n", sep=''))
      res
    }, error = function(e){ cat("ERROR :",conditionMessage(e), "\n") })
  })
  sink()
  
  sink(file.path(o, paste('depths_by_time_chisq_collapsed.txt')))
  r = lapply(sort(unique(df$time))[-1], function(s) {
    tryCatch({
      res = chisq.gof(df %>% mutate(series = Distribution) %>%
                        dplyr::filter(time==s), 'depth', collapse = TRUE)
      cat(paste("\n(t=", s/60, " min. results)\n\n\n", sep=''))
      res
    }, error = function(e){ cat("ERROR :",conditionMessage(e), "\n") })
  })
  sink()
  
  save(df, r, r.imse, file = file.path(o, paste('depths_by_time.RData')))
  
  
  #
  # observed stage duration distributions
  #
  
  df = rbind(
    postpred.samples.obs %>%
      dplyr::filter(max.depth.obs >= deep_dive_depth,
                    depth.start.obs == 1) %>%
      mutate(series = 'Post. Predictive',
             total = length(unique(sample)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      pivot_longer(cols = starts_with("dur."),
                   names_to = 'stage', values_to = 'stage.duration') %>%
      dplyr::filter(!is.na(stage.duration)) %>%
      mutate(stage = fct_recode(stage,
                                'Stage 1' = 'dur.s1',
                                'Stage 2' = 'dur.s2',
                                'Stage 3' = 'dur.s3')) %>%
      group_by(stage, stage.duration, series) %>%
      summarise(prob = n() / total[1],
                n = n(),
                eps = eps[1]) %>%
      group_by(stage, series) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1)),
    validation.obs %>%
      mutate(series = 'Empirical Validation',
             total = length(unique(dive.id)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      pivot_longer(cols = starts_with("dur."),
                   names_to = 'stage', values_to = 'stage.duration') %>%
      mutate(stage = fct_recode(stage,
                                'Stage 1' = 'dur.s1',
                                'Stage 2' = 'dur.s2',
                                'Stage 3' = 'dur.s3')) %>%
      group_by(stage, stage.duration, series) %>%
      summarise(prob = n() / total[1],
                n = n(),
                eps = eps[1]) %>%
      group_by(stage, series) %>%
      mutate(cdf = cumsum(prob),
             cdf.lwr = pmax(cdf - eps, 0),
             cdf.upr = pmin(cdf + eps, 1))
  )
  
  pl = ggplot(df, aes(x = stage.duration/60, y = cdf, ymin = cdf.lwr,
                      ymax = cdf.upr, col = series, fill = series)) +
    geom_ribbon(alpha = .05, col = NA) +
    geom_point() +
    geom_line(lty = 3, alpha = .6) +
    scale_color_brewer('Distribution', type = 'qual', palette = 'Dark2') +
    scale_fill_brewer('Distribution', type = 'qual', palette = 'Dark2') +
    xlab('Observed stage duration (min)') +
    facet_wrap(~stage) +
    ylab('CDF') +
    theme_few()
  
  ggsave(pl, filename = file.path(o, 'stage_duration_cdfs.png'),
         dpi = 'print', width = 16, height = 8)
  
  sink(file.path(o, paste('stage_duration_imse.txt')))
  r.imse = lapply(levels(df$stage), function(s) {
    df2 = df %>% dplyr::filter(stage==s) %>% group_by(series) %>%
      mutate(prob = prob/sum(prob)) %>% ungroup()
    res = imse.gof(df2, 'stage.duration')
    cat(paste("\n(", s," results)\n\n\n", sep=''))
    c(s = s, imse = res)
  })
  sink()
  
  sink(file.path(o, paste('stage_duration_chisq.txt')))
  r = lapply(levels(df$stage), function(s) {
    tryCatch({
      df2 = df %>% dplyr::filter(stage==s) %>% group_by(series) %>%
        mutate(prob = prob/sum(prob)) %>% ungroup()
      res = chisq.gof(df2, 'stage.duration', collapse = FALSE)
      cat(paste("\n(", s," results)\n\n\n", sep=''))
      res
    }, error = function(e){ cat("ERROR :",conditionMessage(e), "\n") })
  })
  sink()
  
  sink(file.path(o, paste('stage_duration_chisq_collapsed.txt')))
  r = lapply(levels(df$stage), function(s) {
    tryCatch({
      df2 = df %>% dplyr::filter(stage==s) %>% group_by(series) %>%
        mutate(prob = prob/sum(prob)) %>% ungroup()
      res = chisq.gof(df2, 'stage.duration', collapse = TRUE)
      cat(paste("\n(", s," results)\n\n\n", sep=''))
      res
    }, error = function(e){ cat("ERROR :",conditionMessage(e), "\n") })
  })
  sink()
  
  sink(file.path(o, paste('stage_duration_ks.txt')))
  r = lapply(levels(df$stage), function(s) {
    tryCatch({
      df2 = df %>% dplyr::filter(stage==s) %>% group_by(series) %>%
        mutate(prob = prob/sum(prob)) %>% ungroup()
      res = ks.gof(df2, 'stage.duration')
      cat(paste("\n(", s," results)\n\n\n", sep=''))
      res
    }, error = function(e){ cat("ERROR :",conditionMessage(e), "\n") })
  })
  sink()
  
  save(df, r, r.imse, file = file.path(o, paste('stage_duration.RData')))
  
  df2 = rbind(
    postpred.samples.obs %>%
      dplyr::filter(max.depth.obs >= deep_dive_depth,
                    depth.start.obs == 1) %>%
      mutate(series = 'Post. Predictive',
             total = length(unique(sample)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      pivot_longer(cols = starts_with(c("dur.s1", 'dur.s2', 'duration.obs')),
                   names_to = 'stage', values_to = 'stage.duration') %>%
      dplyr::filter(!is.na(stage.duration)) %>%
      mutate(stage = fct_recode(stage,
                                'Stage 1' = 'dur.s1',
                                'Stage 2' = 'dur.s2',
                                'Stage 3' = 'dur.s3',
                                'Overall dive' = 'duration.obs')) %>%
      select(stage, series, stage.duration),
    validation.obs %>%
      mutate(series = 'Empirical Validation',
             total = length(unique(dive.id)),
             eps = sqrt(log(2/.05)/(2*total))) %>%
      pivot_longer(cols = starts_with(c("dur.s1", 'dur.s2', 'duration.obs')),
                   names_to = 'stage', values_to = 'stage.duration') %>%
      mutate(stage = fct_recode(stage,
                                'Stage 1' = 'dur.s1',
                                'Stage 2' = 'dur.s2',
                                'Stage 3' = 'dur.s3',
                                'Overall dive' = 'duration.obs')) %>%
      select(stage, series, stage.duration)
  )
  
  pl2 = ggplot(df2, aes(x = series, y = stage.duration/60, col = series)) + 
    facet_wrap(~stage) + 
    geom_boxplot() + 
    ylab('Duration (min.)') + 
    scale_color_brewer('Distribution', type = 'qual', palette = 'Dark2') + 
    theme_few() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank())
  
  ggsave(pl2, filename = file.path(o, 'stage_duration_boxplots.png'),
         dpi = 'print')
  
  ggsave(pl2, filename = file.path(o, paste('stage_duration_boxplots_', 
                                            tag_sex$deployid[tag_id], '.png',
                                            sep = '')),
         dpi = 'print', width = 8, height = 3)
  
  ggsave(pl.pub, 
         filename = file.path(o, paste('depths_by_time_cdf_pub_', 
                                       tag_sex$deployid[tag_id], '.png', 
                                       sep = '')),
         dpi = 'print', width = 8, height = 6)
  
  NULL
}