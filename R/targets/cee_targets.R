cee_targets = list(
  
  tar_target(
    cee_dive_response_targets, 
    c("ZcTag093", "ZcTag095", "ZcTag096", "ZcTag097")
  ),
  
  tar_target(
    cee_surface_response_targets, 
    c("ZcTag069", "ZcTag085")
  ),
  
  tar_target(
    name = cee_dive_response_probs, 
    command = {
      
      # location of MCMC files (also used as output path)
      path = nim_fit
      # path = 'output/mcmc/nim_fit/'
      
      # posterior parameter sample files
      param_sample_files = dir(
        path = path, pattern = 'parameter_samples_[0-9]+', full.names = TRUE
      )
      
      # posterior stage sample files
      stage_sample_files = dir(
        path = path, pattern = 'stage_samples_[0-9]+', full.names = TRUE
      )
      
      # column labels for posterior parameter samples
      param_label_file = dir(
        path = path, pattern = 'parameter_samples_column', full.names = TRUE
      )
      
      # load posterior parameter samples 
      param_samples = do.call(rbind, lapply(param_sample_files, function(f) {
        readRDS(f)
      }))
      
      # load posterior stage samples
      stage_samples = do.call(rbind, lapply(stage_sample_files, function(f) {
        readRDS(f)
      }))
      
      # load data package
      nim_pkg = readRDS(
        file = dir(path = path, pattern = 'nim_pkg', full.names = TRUE)
      )
      
      # label posterior parameter samples
      colnames(param_samples) = readRDS(param_label_file)
      
      burn = 1:2e3
      
      # indices of posterior samples to process
      post_samples = (1:nrow(param_samples))[-burn]
      
      # match tag id in model
      subj_id = which(
        cee_dive_response_targets == nim_pkg$consts$subject_id_labels
      )
      
      # extract pre-exposure data index
      pre_exposure_ind = data.frame(nim_pkg$consts$segments) %>% 
        filter(subject_id == subj_id) %>% 
        select(end_ind) %>% 
        unlist() %>% 
        max()
      
      # subset posterior samples according to parallelization task
      post_samples = post_samples
      
      # draw posterior predictive samples of time to next deep observation
      post_pred_samples = sapply(post_samples, function(sample_ind) {
          
          # extract stage transition coefficients
          betas_tx = array(dim = c(nim_pkg$consts$n_covariates,
                                   nim_pkg$consts$n_stages,
                                   nim_pkg$consts$n_stages - 1))
          for(i in 1:nrow(betas_tx)) {
            for(j in 1:ncol(betas_tx)) {
              for(k in 1:dim(betas_tx)[3]) {
                betas_tx[i,j,k] = param_samples[
                  sample_ind, paste('beta_tx[', i, ', ', j, ', ', k, ']', 
                                    sep = '')
                ]
              }
            }
          }
          
          # forward-simulate dive from predictive distribution
          init_inds = seq(to = pre_exposure_ind, by = 1, length.out = 27)
          fwd_sim = fwd_sim_to_depth(
            stages = stage_samples[sample_ind, init_inds], 
            depths = nim_pkg$data$depths[init_inds], 
            covariates = nim_pkg$data$covariates[, init_inds], 
            n_max = 1e3, 
            nim_pkg = nim_pkg,
            lambda = param_samples[sample_ind, c('lambda[1]', 'lambda[2]')], 
            betas_tx = betas_tx, 
            template_bins = template_bins, 
            times = as.POSIXct(
              x = nim_pkg$data$times[init_inds],
              tz = 'UTC', 
              origin = '1970-01-01 00:00.00 UTC'
            ), 
            timestep = sattag_timestep, 
            lon = cape_hatteras_loc['lon'], lat = cape_hatteras_loc['lat'], 
            depth_threshold = deep_dive_depth
          )
          
          # sampled time to deep depth  
          length(fwd_sim$depths) - 1
          
        })
      
      list(list(
        tag = cee_dive_response_targets,
        samples = post_pred_samples
      ))
      
    }, 
    pattern = map(cee_dive_response_targets), 
    deployment = 'worker',
    memory = 'transient',
    storage = 'worker'
  ),
  
  tar_target(
    name = cee_surface_response_probs, 
    command = {
      
      # location of MCMC files (also used as output path)
      path = nim_fit
      # path = 'output/mcmc/nim_fit/'
      
      # posterior parameter sample files
      param_sample_files = dir(
        path = path, pattern = 'parameter_samples_[0-9]+', full.names = TRUE
      )
      
      # posterior stage sample files
      stage_sample_files = dir(
        path = path, pattern = 'stage_samples_[0-9]+', full.names = TRUE
      )
      
      # column labels for posterior parameter samples
      param_label_file = dir(
        path = path, pattern = 'parameter_samples_column', full.names = TRUE
      )
      
      # load posterior parameter samples 
      param_samples = do.call(rbind, lapply(param_sample_files, function(f) {
        readRDS(f)
      }))
      
      # load posterior stage samples
      stage_samples = do.call(rbind, lapply(stage_sample_files, function(f) {
        readRDS(f)
      }))
      
      # load data package
      nim_pkg = readRDS(
        file = dir(path = path, pattern = 'nim_pkg', full.names = TRUE)
      )
      
      # label posterior parameter samples
      colnames(param_samples) = readRDS(param_label_file)
      
      burn = 1:2e3
      
      # indices of posterior samples to process
      post_samples = (1:nrow(param_samples))[-burn]
      
      # match tag id in model
      subj_id = which(
        cee_surface_response_targets == nim_pkg$consts$subject_id_labels
      )
      
      # extract pre-exposure data index
      pre_exposure_ind = data.frame(nim_pkg$consts$segments) %>% 
        filter(subject_id == subj_id) %>% 
        select(end_ind) %>% 
        unlist() %>% 
        max()
      
      # subset posterior samples according to parallelization task
      post_samples = post_samples
      
      # match tag id in collection of raw sattags
      raw_sattag_ind = which(
        cee_surface_response_targets == sapply(raw_sattags, function(x) x$tag)
      )
      
      # extract matching raw sattag and id of exposed dive
      raw_tag = raw_sattags[[raw_sattag_ind]]
      exposed_dive_id = raw_tag$diveIds[min(which(raw_tag$exposed == 1))]
      
      # last depth that can be associated with exposed dive
      last_exposed_depth = raw_tag$depths[
        max(which(raw_tag$diveIds == exposed_dive_id)) + 1
      ]
      
      # draw posterior predictive samples of time to next deep observation
      post_pred_samples = sapply(post_samples, function(sample_ind) {
        
        # extract stage transition coefficients
        betas_tx = array(dim = c(nim_pkg$consts$n_covariates,
                                 nim_pkg$consts$n_stages,
                                 nim_pkg$consts$n_stages - 1))
        for(i in 1:nrow(betas_tx)) {
          for(j in 1:ncol(betas_tx)) {
            for(k in 1:dim(betas_tx)[3]) {
              betas_tx[i,j,k] = param_samples[
                sample_ind, paste('beta_tx[', i, ', ', j, ', ', k, ']', 
                                  sep = '')
              ]
            }
          }
        }
        
        # forward-simulate dive from predictive distribution
        init_inds = seq(to = pre_exposure_ind, by = 1, length.out = 27)
        fwd_sim = fwd_sim_to_depth(
          stages = stage_samples[sample_ind, init_inds], 
          depths = nim_pkg$data$depths[init_inds], 
          covariates = nim_pkg$data$covariates[, init_inds], 
          n_max = 1e3, 
          nim_pkg = nim_pkg,
          lambda = param_samples[sample_ind, c('lambda[1]', 'lambda[2]')], 
          betas_tx = betas_tx, 
          template_bins = template_bins, 
          times = as.POSIXct(
            x = nim_pkg$data$times[init_inds],
            tz = 'UTC', 
            origin = '1970-01-01 00:00.00 UTC'
          ), 
          timestep = sattag_timestep, 
          lon = cape_hatteras_loc['lon'], lat = cape_hatteras_loc['lat'], 
          depth_threshold = last_exposed_depth,
          deeper = FALSE
        )
        
        # sampled time to last depth associated with dive
        length(fwd_sim$depths) - 1
        
      })
      
      list(list(
        tag = cee_surface_response_targets,
        samples = post_pred_samples
      ))
      
    }, 
    pattern = map(cee_surface_response_targets), 
    deployment = 'worker',
    memory = 'transient',
    storage = 'worker'
  ),
  
  tar_target(
    name = cee_dive_response_summaries, 
    command = {
    
      depth_threshold = deep_dive_depth
      
      # extract data for exposure conditions
      exposure_contexts = lapply(raw_sattags, function(tag) {
        
        # skip processing if tag was never exposed
        if(!is.finite(tag$exposure_time)) {
          return(NULL)
        }
        
        # last observation before exposure
        pre_exposed_ind = min(which(tag$exposed == 1)) - 1
        
        # 
        # extract dive context immediately before exposure
        #
        
        # proportion of recent observations at depth, as a covariate
        prop_recent_deep = {
          # window at which recent observations begins
          window_start = tag$times[pre_exposed_ind] - 
            duration(135, units = 'minutes')
          # data indices of recent observations
          past_inds = 1:(pre_exposed_ind-1)
          window_inds = past_inds[window_start <= tag$times[past_inds]]
          # compute proportion of recent observations spent below a depth
          sum(tag$depths[window_inds] >= depth_threshold) / length(window_inds)
        }
        
        #
        # extract time until next deep depth observation
        #
        
        
        # data index at which first deep depth post-exposure is observed
        next_deep_obs = which(tag$exposed == 1)[
          min(which(tag$depths[tag$exposed == 1] >= depth_threshold))
        ]
        
        # censored data due to gaps in data
        if(any(tag$gap_after[(pre_exposed_ind + 1):next_deep_obs])) {
          next_deep_obs = Inf
        }
        
        #
        # extract time until the last depth associated with exposed dive
        #
        
        # extract id of exposed dive
        exposed_dive_id = tag$diveIds[min(which(tag$exposed == 1))]
        
        # index of last depth that can be associated with exposed dive
        last_exposed_depth_ind = max(which(tag$diveIds == exposed_dive_id)) + 1
        
        list(
          tag = tag$tag,
          depth.bin.prev = tag$depth.bin[pre_exposed_ind-1],
          depth.bin = tag$depth.bin[pre_exposed_ind],
          depth = tag$depths[pre_exposed_ind],
          prop_recent_deep = prop_recent_deep,
          daytime = tag$daytime[pre_exposed_ind],
          moonlit = tag$moonlit[pre_exposed_ind],
          time_to_deep = next_deep_obs - pre_exposed_ind,
          time_to_end = last_exposed_depth_ind - pre_exposed_ind,
          pre_exposure_time = tag$times[pre_exposed_ind]
        )
      })
      
      
      #
      # look at dive response
      #
      
      df = do.call(rbind, lapply(cee_dive_response_probs, 
                                 function(post_samples) {
        
        # find exposure context associated with posterior samples
        tag_id = which(
          post_samples$tag == sapply(exposure_contexts, function(x) x$tag)
        )
        
        # compute p-value (post. pred. prob. of less time to deep depth)
        p_val = mean(
          post_samples$samples <= exposure_contexts[[tag_id]]$time_to_deep
        )
        
        # package results
        data.frame(
          tag = post_samples$tag, 
          exposure_depth_bin = exposure_contexts[[tag_id]]$depth.bin,
          daytime = exposure_contexts[[tag_id]]$daytime,
          moonlit = exposure_contexts[[tag_id]]$moonlit,
          prop_recent_deep = exposure_contexts[[tag_id]]$prop_recent_deep,
          time_to_deep = exposure_contexts[[tag_id]]$time_to_deep,
          p = round(p_val,2)
        )
      }))
      
      # package results
      list(
        exposure_contexts = exposure_contexts,
        pvals = df
      )
      
    }
  ),
  
  tar_target(
    name = exposure_context_plots,
    command = {

      zc93_ind = which(sapply(raw_sattags, function(res) res$tag) == 'ZcTag093')
      zc95_ind = which(sapply(raw_sattags, function(res) res$tag) == 'ZcTag095')
      
      
      plot_tag_segment = function(tag, pre_len, post_len) {

        # determine when tag was exposed
        exposed_ind = min(which(tag$exposed == 1))
        
        # times before and after exposure to study
        window_start = tag$times[exposed_ind] - 
          duration(pre_len, units =  'hours')
        window_end = tag$times[exposed_ind] + 
          duration(post_len, units =  'hours')
        
        # data associated with window indices
        window_inds = which(
          (tag$times >= window_start) & (tag$times <= window_end)
        )
        
        # template bins depth ranges
        bin_ranges = cbind(
          lwr = template_bins$center - template_bins$halfwidth,
          upr = template_bins$center + template_bins$halfwidth
        )
        
        df = data.frame(
          time = tag$times[window_inds],
          depth = tag$depths[window_inds],
          depth_lwr = bin_ranges[tag$depth.bin[window_inds], 'lwr'],
          depth_upr = bin_ranges[tag$depth.bin[window_inds], 'upr']
        )
        
        ggplot(df, 
               aes(x = time, y = depth, ymin = depth_lwr, ymax = depth_upr)) + 
          geom_pointrange() + 
          geom_line() + 
          geom_vline(xintercept = tag$exposure_time, lty = 3) + 
          geom_hline(yintercept = deep_dive_depth, lty = 3) + 
          scale_y_reverse('Depth (m)', limits = rev(c(0, 1500))) + 
          theme_few()
      }
      
      pre_len = 2.5
      
      pl93 = plot_tag_segment(
        tag = raw_sattags[[zc93_ind]], pre_len = pre_len, post_len = 1
      ) + ggtitle('Zc093') + xlab('')
      
      pl95 = plot_tag_segment(
        tag = raw_sattags[[zc95_ind]], pre_len = pre_len, post_len = 1
      ) + ggtitle('Zc095') + xlab('Time (UTC)')
      
      # save file
      f = file.path('output', 'cee')
      dir.create(f, showWarnings = FALSE, recursive = TRUE)
      ggsave(
        ggarrange(pl93, pl95, nrow = 2), 
        filename = file.path(f, paste(tar_name(), '.png', sep = '')),
        dpi = 'print', width = 8, height = 6
      )
      
    }
  )
  
)
