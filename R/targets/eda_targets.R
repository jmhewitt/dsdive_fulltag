eda_targets = list(
  
  tar_target(
    name = mle_speeds, 
    command = {
      
      # unwrap tag
      tag = raw_sattags[[1]]
      
      # depth bin widths and count
      n_bins = nrow(template_bins)
      bin_widths = 2 * template_bins$halfwidth
      
      # 
      nobs = 1e2
      
      tick = proc.time()[3]
      mle_params = do.call(rbind, lapply(1:nobs, function(ind) {
        
        # time between observations (sec)
        tstep = diff(as.numeric(tag$times[ind + 0:1]))
        
        # depth bin at start and end of observation
        s0 = tag$depth.bin[ind]
        sf = tag$depth.bin[ind + 1]
        
        # # optimize transition parameters
        # o = optim(c(0,0), fn = function(theta) {
        #   expm(tstep * buildInfinitesimalGenerator(
        #     pi = plogis(theta[1]), 
        #     lambda = exp(theta[2]), 
        #     M = n_bins, 
        #     stage = 6, # unrestricted bin transitions
        #     widths = bin_widths
        #   ))[s0, sf]
        # }, control = list(fnscale = -1, maxit = 1e3), method = 'BFGS')
        
        piseq = c(.01, .5, .99)
        
        # profile-optimization of transition parameters
        o = lapply(piseq, function(pi) {
          optim(0,  fn = function(log_lambda) {
            expm(tstep * buildInfinitesimalGenerator(
              pi = pi, 
              lambda = exp(log_lambda), 
              M = n_bins, 
              stage = 6, # unrestricted bin transitions
              widths = bin_widths
            ))[s0, sf]
          }, control = list(fnscale = -1), method = 'BFGS')
        })
        
        if(any(sapply(o, function(o) o$convergence) != 0)) {
          warning(paste('Convergence failed for obs_ind', ind))
        }
        
        o_ind = which.max(sapply(o, function(o) o$value))
        
        # extract and back-transform parameters
        data.frame(pi = piseq[o_ind], lambda = exp(o[[o_ind]]$par))
      }))
      tock = proc.time()[3]
      
      tock-tick
      
      mle_params$t = tag$times[1:nobs]
      mle_params$ind = 1:nobs
      mle_params$depth = tag$depths[1:nobs]
      mle_params$depth.next = tag$depths[1:nobs + 1]
      mle_params$apparent_speed = abs(mle_params$depth.next - mle_params$depth)/300
      mle_params$depth.bin = tag$depth.bin[1:nobs]
      mle_params$depth.bin.next = tag$depth.bin[1:nobs + 1]
      mle_params$depth.bin.prev = c(NA, tag$depth.bin[1:(nobs-1)])
      mle_params$depth.bin.txd = abs(mle_params$depth.bin.next - mle_params$depth.bin)
        
      ggplot(mle_params %>% pivot_longer(cols = c(pi,lambda, depth), 
                                         names_to = 'param'), 
             aes(x = t, y = value)) + 
        geom_line() + 
        geom_point() + 
        theme_few() + 
        facet_wrap(~param, scales = 'free', ncol = 1)
      
      
      library(viridis)
      
      ggplot(mle_params %>% 
               mutate(pi = factor(pi),
                      lambda = cut(lambda, breaks = c(-.1,0.5,1,6))) %>%
               pivot_longer(cols = c(pi,lambda), names_to = 'parameter'), 
             aes(x = t, y = depth, 
                 group = 1,
                 col = value)) +
                 # col = factor(pi))) +
                 # col = cut(lambda, breaks = 5))) + 
        facet_wrap(~parameter, ncol = 1) + 
        geom_line() + 
        geom_point() + 
        scale_y_reverse() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2') +
        theme_few() 
      
      
      ggplot(mle_params, aes(x = depth.bin.next, y = lambda)) + 
        geom_point() + 
        scale_x_continuous(breaks = 1:16) 
      
      ggplot(mle_params, aes(x = factor(depth.bin), y = lambda,
                             col = depth.bin >= depth.bin.prev)) + 
        geom_boxplot() +
        geom_point() 
      
      ggplot(mle_params, aes(x = factor(pi), y = lambda)) + 
        geom_boxplot() +
        geom_point() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2')
      
      ggplot(mle_params, aes(x = factor(pi), y = lambda),
                             col = depth.bin >= depth.bin.prev) + 
        geom_boxplot() +
        geom_point() + 
        scale_color_brewer(type = 'qual', palette = 'Dark2')
      
      
      browser()
    }, 
    pattern = map(raw_sattags)
  ),
  
  tar_target(
    name = recovery_eda_plot,
    command = {
      
      tar_load(nim_pkg)
      tar_load(template_bins)
      tar_load(deep_dive_depth)
      
      # forward-looking time windows in which to look for deep depths
      hrs = 1:6
      
      # extract data to explore trends in how covariates are associated with 
      # future observations of deep depths, process by segment
      df = do.call(rbind, lapply(1:nim_pkg$consts$n_segments, function(seg_id) {
        
        # data-indices associated with segment
        seg_inds = seq(
          from = nim_pkg$consts$segments[seg_id, 'start_ind'],
          to = nim_pkg$consts$segments[seg_id, 'end_ind'], 
          by = 1
        )
        
        # # exploratory data, by observation within segment
        # do.call(rbind, lapply(seg_inds, function(ind) {
        #   data.frame(
        #     # depth of current observation
        #     start_depth = nim_pkg$data$depths[ind],
        #     # covariates for current observation
        #     t(nim_pkg$data$covariates[,ind]),
        #     # copy-paste the exploratory time windows being studied
        #     hr = hrs,
        #     # identify if a deep depth is observed in the near-future from ind
        #     any_deep = sapply(hrs, function(hr) {
        #       # ideal data indices to investigate
        #       nominal_inds = ind + 1:(12*hr)
        #       # restrict indices to lie entirely within the segment
        #       tgt_inds = intersect(seg_inds, nominal_inds)
        #       # skip processing if we are near the end of the segment
        #       if(length(nominal_inds) != length(tgt_inds)) {
        #         NA
        #       }
        #       # check to see if there are deep observations in the near-future
        #       any(
        #         template_bins$center[nim_pkg$data$depths[tgt_inds]] > 
        #           deep_dive_depth
        #       )
        #     })
        #   )
        # }))
        
        # exploratory data, by observation within segment
        do.call(rbind, lapply(seg_inds, function(ind) {
          # identify if a deep depth is observed in the near-future from ind
          do.call(rbind, lapply(hrs, function(hr) {
            # ideal data indices to investigate
            nominal_inds = ind + 1:(12*hr)
            # restrict indices to lie entirely within the segment
            tgt_inds = intersect(seg_inds, nominal_inds)
            # "skip" processing if we are near the end of the segment
            if(length(nominal_inds) != length(tgt_inds)) {
              any_deep = NA
              tmp_covariates = NA * t(nim_pkg$data$covariates[,ind])
            } else {
              any_deep = any(
                template_bins$center[nim_pkg$data$depths[tgt_inds]] > 
                  deep_dive_depth
              )
              tmp_covariates = t(nim_pkg$data$covariates[,tail(tgt_inds, 1)])
            }
            # build results
            data.frame(
              # depth of first observation in window
              start_depth = nim_pkg$data$depths[ind],
              # covariates for first observation in window
              t(nim_pkg$data$covariates[,ind]),
              # exploratory time window being studied
              hr = hr,
              # check to see if there are deep observations in the near-future
              any_deep = any_deep,
              # covariates for final observation in target window
              tmp_covariates
            )
          }))
        }))
        
      }))
      
      # empirical estimates for P(deep_depth | data)
      df.eda = df %>% 
        filter(
          # only analyze EDA windows with valid data
          is.finite(any_deep), 
          # restrict EDA to observations that start and end during the day
          daytime == 1, 
          daytime.1 == 1,
          # restrict EDA to windows that start on surface
          start_depth == 1,
          # analyze shorter prediction windows, which have non-trivial dist'n.s
          hr <= 3) %>%
        # convert bool's to indicators
        mutate(any_deep = as.numeric(any_deep)) %>% 
        # group data to allow for empirical proportion estimates
        group_by(prop_recent_deep, hr) %>% 
        # compute empirical proportions, and extract data for naive CI's
        summarise(prob_any_deep = mean(any_deep),
                  x = sum(any_deep),
                  n = n()) %>% 
        # simplify data structure
        ungroup()
      
      # enrich empirical estimates for P(deep_depth | data) with naive CIs
      df.eda = cbind(
        df.eda, 
        do.call(rbind, apply(df.eda, 1, function(r) {
          res = prop.test(x = r['x'], n = r['n'])
          data.frame(prop_lwr = res$conf.int[1], prop_upr = res$conf.int[2])
        }))
      )
      
      pl = ggplot(df.eda %>% 
                    mutate(hr = paste(hr, 'h window', sep = '')), 
                  aes(x = prop_recent_deep, y = prob_any_deep, 
                      ymin = prop_lwr, ymax = prop_upr)) + 
        geom_pointrange() + 
        facet_wrap(~hr) + 
        ylab(expression(hat(p))) + 
        xlab('Proportion of previous hour observed at deep depths (>800m)') + 
        # ggtitle('Empirical prob. to observe deep depth in next window') + 
        theme_few() +
        theme(axis.title.y = element_text(angle = 0, vjust = .5))
      
      # save file
      f = file.path('output', 'eda')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      ggsave(pl, filename = file.path(f, paste(tar_name(), '.png')), 
             dpi = 'print')
      
      # eda is fine, but it is not able to account for the effects of the 
      # movement type, so we model this especially because this is useful for 
      # conditioning on the movement type that a whale was doing at time of 
      # exposure.  the modeling also lets us borrow strength from observations
      # made at depth bin 2, etc., so we get a more broad use of the data.
      # 
      # also i think that the most important part of the study is to see that 
      # the numbers fit roughly well during the first hour post-exposure.  this 
      # should be something that we can check using out-of-sample validation 
      # without too much difficulty
      #
      # more differences is that it is harder to further subset the data s.t.
      # we only analyze daytime observations... this ends up throwing away a 
      # lot of data that we could have otherwise modeled via methods that 
      # borrow strength
      
      
    }
  ),
  
  
  tar_target(
    name = dive_survival_eda_plot,
    command = {
      
      # subset validation dataset
      df.eda.raw = validation_df_deep_surv %>% 
        # allow indexing back to the order of the validation dives
        mutate(val_ind = 1:n()) 
      
      # discrete event times to study
      tseq = seq(
        from = 0,
        to = max(pmin(df.eda.raw$deep_time, df.eda.raw$nobs)),
        by = 1
      )
      
      # simplify conditioning covariate
      df.eda.raw$prop_recent_deep = cut(
        x = df.eda.raw$prop_recent_deep,
        breaks = (0:8)/8,
        include.lowest = TRUE
      )
      
      # empirical estimates of survival probabilities
      df.eda = do.call(rbind, lapply(tseq, function(event_time) {
        df.eda.raw %>% 
          # portion of dataset with sufficiently long observation periods
          filter(nobs >= event_time) %>% 
          # group by remaining covariate of interest
          group_by(prop_recent_deep) %>%
          # empirical probability estimate
          summarise(
            # number of observations and events in data
            nobs = n(),
            nevents = sum(deep_time <= event_time),
            # proportion of deep dives by event_time
            p = nevents / nobs,
            # record event time, for plotting
            event_time = event_time
          )
      }))
      
      # enrich empirical estimates  with naive CIs
      df.eda = cbind(
        df.eda, 
        do.call(rbind, apply(df.eda, 1, function(r) {
          res = prop.test(x = as.numeric(r['nevents']), 
                          n = as.numeric(r['nobs']))
          data.frame(p_lwr = res$conf.int[1], p_upr = res$conf.int[2])
        }))
      )
      
      merge_tables = function(x1, x2) {
        # merge tables to get complete list of index values
        m = x1 %>% full_join(x2, by = 'x')
        # set na's in frequency counts to 0
        m[is.na(m)] = 0
        # return aggregated counts
        m %>%
          mutate(Freq = Freq.x + Freq.y) %>% 
          select(x, Freq)
      }
      
      # combine posterior predictive cdfs by prop_recent_deep
      df.post = do.call(
        rbind, lapply(levels(df.eda$prop_recent_deep), function(prop) {
          # ids of validation dives with prop covariate
          val_ids = unique(df.eda.raw$val_ind[
            df.eda.raw$prop_recent_deep == prop
          ])
          # aggregate pmfs as a mixture
          pmf = validate_deep_surv_distn[[val_ids[1]]]
          if(length(val_ids) > 1) {
            for(vid in val_ids[-1]) {
              pmf = merge_tables(
                x1 = pmf, x2 = validate_deep_surv_distn[[vid]]
              )
            }
            # normalize distribution assuming all component pmfs are equally wtd.
            pmf$Freq = pmf$Freq / length(val_ids)
            # sort pmf and compute cdf
            pmf = pmf[order(pmf$x),]
            pmf$cdf = cumsum(pmf$Freq)
            # return result
            data.frame(prop_recent_deep = prop, pmf)
          }
      }))
      
      # munge names and update factor order
      df.post = df.post %>% mutate(event_time = x, p = cdf)
      df.post$prop_recent_deep = factor(
        x = df.post$prop_recent_deep, 
        levels = levels(df.eda$prop_recent_deep)
      )
      
      pl = ggplot(df.eda, 
                  aes(y = p, x = event_time)) + 
        # validation distribution
        geom_line() +
        geom_pointrange(aes(ymin = p_lwr, ymax = p_upr)) +
        # posterior predictive distribution
        geom_line(data = df.post %>% filter(event_time %in% tseq),
                  col = 'salmon') +
        # formatting
        facet_wrap(~prop_recent_deep,) +
        ylab('P(Deep depth)') + 
        xlab('Num. observations forward') + 
        theme_few() +
        theme(panel.grid.major = element_line(colour = 'grey95'))
      
      # save file
      f = file.path('output', 'eda')
      dir.create(path = f, showWarnings = FALSE, recursive = TRUE)
      ggsave(pl, filename = file.path(f, paste(tar_name(), '.png')), 
             dpi = 'print', width = 12, height = 8)
      
      # This is a naive estimator of the survivor function, which can be 
      # improved by using a Kaplan-Meier curve instead, which accounts for 
      # censoring.  however, the improvement is perhaps not important because 
      # we still have the issue that the underlying data are highly dependent, 
      # as we take multiple sequential observations from a single individual,
      # so a proper survival analysis would need to take this into account for 
      # estimation.  the naive estimator is nice, in particular, because it 
      # highlights how the data is partitioned.
      #
      # a process-based estimator does a good job of accounting for this 
      # aspect of the data, while also allowing for the possibility to conduct 
      # other analyses (not just survival analyses), and it has a generative
      # structure that makes it easier to intuitively check the model structure,
      # whereas survival analysis may be a little less direct.
      
    }
  )
  
)
