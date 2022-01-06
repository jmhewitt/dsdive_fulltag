segment_dives = function(dive_label_plot_dir, label_diagnostic_plot_dir, 
                         template_bins, exploratory_merge_ratios, 
                         depth_files, merge_ratio, tag_info, deep_threshold,
                         timestep) {
  
  #
  # find and process raw data
  #
  
  # process raw files
  tag.records = lapply(depth_files, function(f) {
    
    # load data
    d = read.csv(file = f)
    d$Date = as.POSIXct(d$Date, tz = 'UTC', origin = '1970-01-01 00:00.00 UTC')
    
    # extract tag name
    tag.name = str_extract(f, 'Zc[0-9A-Za-z]+')
    
    # map all depths to standardized bins
    d$depth.bin = sapply(d$Depth, function(depth) {
      which.min(abs(depth - template_bins$center))
    })
    d$depth.standardized = template_bins$center[d$depth.bin]
    
    # label dives using range of ratios
    label.seq = lapply(exploratory_merge_ratios, function(mr) {
      list(
        labels = dive.segmentation(y = d$depth.bin, merge.ratio = mr,
                                   depth.bins = template_bins, times = d$Date,
                                   timestep = timestep),
        merge.ratio = mr
      )
    })
    
    # package results
    list(
      depths = d,
      label.seq = label.seq,
      name = tag.name
    )
  })
  
  
  #
  # segmentation diagnostics
  #
  
  # extract summary features from records
  labels.diagnostics = do.call(rbind, lapply(tag.records, function(record) {
    # loop over segmentations within record
    do.call(rbind, lapply(record$label.seq, function(segmentation) {
      # join depth data with segmentation labels and configuration
      cbind(record$depths, 
            mode = segmentation$labels, 
            merge.ratio = segmentation$merge.ratio) %>% 
        # restrict summary to known diving behavior
        filter(mode > 0) %>%
        # compute dive-level summaries
        group_by(mode, merge.ratio) %>% 
        summarise(maxd = max(depth.standardized),
                  start = min(Date),
                  end = max(Date),
                  duration = difftime(time1 = end, time2 = start, 
                                      units = 'mins'),
                  deep = (maxd >= deep_threshold)) %>% 
        # aggregate summaries by dive type
        group_by(deep) %>% 
        summarise(n = length(deep), 
                  duration.q25 = quantile(duration, probs = .25),
                  duration.q5 = quantile(duration, probs = .5),
                  duration.mean = mean(duration),
                  duration.q75 = quantile(duration, probs = .75),
                  merge.ratio = merge.ratio[1]) %>% 
        # add tag label and format for plotting
        mutate(tag = record$name,
               deep = factor(deep, 
                             levels = c(TRUE, FALSE), 
                             labels = c('Deep',' Shallow'))
        )
    }))
  }))
  
  
  # compare number of dives wrt. merge ratios
  pl = ggplot(labels.diagnostics %>% 
                group_by(merge.ratio, deep) %>%
                summarise(n.mean = mean(n),
                          n.lwr = quantile(n, probs = .25),
                          n.upr = quantile(n, probs = .75)), 
              aes(x = merge.ratio, y = n.mean)) + 
    geom_line() + 
    geom_line(mapping = aes(x = merge.ratio, y = n.lwr), lty = 3, 
              inherit.aes = FALSE) + 
    geom_line(mapping = aes(x = merge.ratio, y = n.upr), lty = 3, 
              inherit.aes = FALSE) + 
    xlab('Merge ratio') + 
    ylab('Dives per record') + 
    facet_wrap(~deep, ncol = 1, scales = 'free') + 
    theme_few() + 
    theme(panel.background = element_blank(), strip.placement = 'inside') + 
    ggtitle('(Solid: mean across tags; Dotted: 25% and 75% quantiles)')
  
  ggsave(pl, filename = file.path(label_diagnostic_plot_dir, 
                                  'dives_per_ratio.pdf'))
  
  
  
  # compare dive durations wrt. merge ratios
  pl = ggplot(labels.diagnostics %>% 
                group_by(merge.ratio, deep) %>% 
                summarise(duration.mean = mean(duration.mean),
                          duration.q25 = quantile(duration.q25, probs = .25),
                          duration.q75 = quantile(duration.q75, probs = .75)), 
              aes(x = merge.ratio, y = duration.mean)) + 
    geom_line() + 
    geom_line(mapping = aes(x= merge.ratio, y = duration.q25), 
              inherit.aes = FALSE, lty = 3) + 
    geom_line(mapping = aes(x= merge.ratio, y = duration.q75), 
              inherit.aes = FALSE, lty = 3) + 
    xlab('Merge ratio') + 
    ylab('Dive duration (min)') + 
    facet_wrap(~deep, ncol = 1, scales = 'free') + 
    theme_few() + 
    theme(panel.background = element_blank(), strip.placement = 'inside') + 
    ggtitle('(Solid: mean across tags; Dotted: Est. 25% and 75% quantiles)')
  
  ggsave(pl, filename = file.path(label_diagnostic_plot_dir, 
                                  'durations_per_ratio.pdf'))
  
  
  #
  # tag plots and save labels
  #
  
  # merge ratios selected
  merge.ratios.ideal = rep(merge_ratio, length(tag.records))
  
  dive_labels = mapply(FUN = function(record, ratio) {
    
    # find the labeled depth series that best matches the desired merge ratio
    merge.ind = which.min(abs(
      sapply(record$label.seq, function(lab) lab$merge.ratio) - ratio
    ))
    
    # extract dive labels
    labs = record$label.seq[[merge.ind]]$labels
    
    # plot dive record
    pl = tagplot(depths = record$depths, depth.bins = template_bins, 
                 dives.labeled = labs, cee_starts = tag_info$cee_start, 
                 depth_mark = deep_threshold)
    
    # save plot of dive record 
    ggsave(pl, filename = file.path(dive_label_plot_dir, 
                                    paste(record$name, '.pdf', sep = '')),
           dpi = 'print', height = 12, width = 12*12, limitsize = FALSE)
    
    
    # return dive labels
    list(list(labels = labs, tag = record$name))
    
  }, tag.records, merge.ratios.ideal)

  dive_labels
}

