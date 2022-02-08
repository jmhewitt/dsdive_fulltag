validation_nstep_report_script = tar_target(
  name = validation_nstep_report, 
  command = {
    
    # combine scores from all validation targets
    scores = do.call(
      rbind, lapply(validation_nstep_processing, function(x) {
        x$scores$from_depth_bin = x$from_depth_bin
        x$scores
      })
    )
    
    # summarize overall prediction error, 1/2-step
    summary(scores)
    
    # error does not seem to be sensitive to amount of data used to est. xf
    ggplot(scores, aes(x = latent_state_training_inds, y = crps_bin_1)) + 
      geom_point() + 
      stat_smooth() + 
      scale_x_log10() + 
      theme_few()
    
    # error tends to decrease as starting depth bin increases
    scores %>% 
      group_by(from_depth_bin) %>% 
      summarise(
        crps_bin_1 = mean(crps_bin_1),
        crps_bin_2 = mean(crps_bin_2)
      ) %>% 
      left_join(
        data.frame(depth_bin = 1:data_pkg$consts$n_bins,
                   width = data_pkg$consts$widths),
         by = c(from_depth_bin = 'depth_bin')
      ) %>% 
      ggplot(aes(x = width, y = crps_bin_1)) + 
      geom_point() + 
      theme_few() +
      stat_smooth(method = 'lm')
    
    #
    # compare empirical transition matrix to modeled transition matrix
    #
    
    # observed depth bin transitions
    empirical = do.call(rbind, lapply(validation_nstep_processing, function(x) {
      data.frame(from_depth_bin = x$from_depth_bin,
                 to_depth_bin_1 = x$to_depth_bin_1,
                 to_depth_bin_2 = x$to_depth_bin_2)
    }))
    
    # empirical distribution of 1-step depth bin transition matrix
    empirical_1_step = do.call(rbind, lapply(1:data_pkg$consts$n_bins, 
                                             function(from_bin) {
      dist = tabulate(
        bin = empirical %>% 
          filter(from_depth_bin == from_bin) %>% 
          select(to_depth_bin_1) %>% 
          unlist(), 
        nbins = data_pkg$consts$n_bins)
      dist / sum(dist)
    }))
    
    # sample size for each row of empirical transition matrix
    empirical_n = do.call(rbind, lapply(1:data_pkg$consts$n_bins, 
                                        function(from_bin) {
      dist = tabulate(
        bin = empirical %>% 
          filter(from_depth_bin == from_bin) %>% 
          select(to_depth_bin_1) %>% 
          unlist(), 
        nbins = data_pkg$consts$n_bins)
     sum(dist)
   }))
    
    # empirical distribution of 2-step depth bin transition matrix
    empirical_2_step = do.call(rbind, lapply(1:data_pkg$consts$n_bins, 
                                             function(from_bin) {
     dist = tabulate(
       bin = empirical %>% 
         filter(from_depth_bin == from_bin) %>% 
         select(to_depth_bin_2) %>% 
         unlist(), 
       nbins = data_pkg$consts$n_bins)
     dist / sum(dist)
   }))
    
  pred_1_step = do.call(rbind, lapply(1:data_pkg$consts$n_bins, 
                                      function(from_bin) {
    colMeans(do.call(rbind, lapply(validation_nstep_processing, 
                                   function(val_pack) {
      if(val_pack$from_depth_bin == from_bin) {
        return(val_pack$pred_depth_bin_1)
      } else {
        return(NULL)
      }
    })))
  }))
  
  pred_2_step = do.call(rbind, lapply(1:data_pkg$consts$n_bins, 
                                      function(from_bin) {
  colMeans(do.call(rbind, lapply(validation_nstep_processing, 
                                 function(val_pack) {
                                   if(val_pack$from_depth_bin == from_bin) {
                                     return(val_pack$pred_depth_bin_2)
                                   } else {
                                     return(NULL)
                                   }
                                 })))
}))
  
  
  txmat_to_plottable = function(x) {
    data.frame(from_depth_bin = 1:data_pkg$consts$n_bins, x) %>% 
      pivot_longer(cols = starts_with('X'), names_to = 'to_depth_bin', 
                   values_to = 'p') %>%
      mutate(to_depth_bin = as.numeric(gsub('X', '', to_depth_bin))) %>% 
      group_by(from_depth_bin) %>%
      arrange(to_depth_bin) %>% 
      mutate(cdf = cumsum(p)) %>%
      ungroup()
  }
  
  rbind(
    cbind(txmat_to_plottable(empirical_1_step), Method = 'Empirical') %>% 
      left_join(
        data.frame(from_depth_bin = 1:data_pkg$consts$n_bins,
                   n = empirical_n),
        by = 'from_depth_bin') %>% 
      mutate(hw = 1.96 * sqrt(p * (1-p) / n),
             lwr = pmax(p - hw, 0), upr = pmin(p + hw, 1)) %>%
      select(-hw, -n)
    ,
    cbind(txmat_to_plottable(pred_1_step), Method = 'Modeled') %>% 
      mutate(lwr = p, upr = p)
  ) %>% 
    ggplot(aes(x = to_depth_bin, y = cdf, col = Method, 
               ymin = lwr, ymax = upr)) + 
    geom_line() + 
    geom_point() +
    # facet_wrap(~from_depth_bin, scales = 'free') +
    facet_wrap(~from_depth_bin) + 
    scale_color_brewer(type = 'qual', palette = 'Dark2') +
    scale_x_continuous(breaks = 1:data_pkg$consts$n_bins) + 
    theme_few()
  
  
  #
  #
  #
  
  # ways to decompose prediction error
  #  - function of "proportion of time" covariate
  
  scores$prop_recent_val = sapply(validation_nstep_processing, function(val_pack) {
    mean(template_bins$center[data_pkg$data$depth_bins[
      seq(to = val_pack$validation_nstep_ind, length.out = 12)
    ]] > 800)
  })
  
  ggplot(scores %>% 
           group_by(prop_recent_val) %>% 
           summarise(crps_bin_1 = mean(crps_bin_1)), 
         aes(x = prop_recent_val, y = crps_bin_1)) + 
    geom_point() + 
    theme_few()
  
  
  
  browser()
  
    
  }
)