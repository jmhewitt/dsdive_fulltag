plot_targets = list(

  
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
      
  }),
  
  tar_target(
    name = depth_ctmc_illustration,
    command = {
      
      zc93_ind = which(sapply(raw_sattags, function(res) res$tag) == 'ZcTag093')
      
      plot_tag_segment = function(tag, pre_len, post_len, mvmt_type) {
        # Parameters: 
        #  mvmt_type - 'fast_descent', 'slow_descent', 'fast_ascent', etc...
        
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
        
        nim_pkg = readRDS(file.path(
          'output', 'mcmc', 'fit_marginalized_model', 'nim_pkg.rds'
        ))
        
        ind = which(rownames(nim_pkg$consts$stage_defs) == mvmt_type)
        
        Pmat = expm::expm(
          nim_pkg$consts$tstep * buildInfinitesimalGenerator(
            pi = nim_pkg$consts$pi[nim_pkg$consts$stage_defs[ind, 1]],
            lambda = nim_pkg$inits$lambda[nim_pkg$consts$stage_defs[ind, 2]],
            M = nim_pkg$consts$n_bins,
            stage = 6,
            widths = nim_pkg$consts$widths
          )
        )
        
        df = data.frame(
          time = tag$times[window_inds],
          depth = tag$depths[window_inds],
          depth_lwr = bin_ranges[tag$depth.bin[window_inds], 'lwr'],
          depth_upr = bin_ranges[tag$depth.bin[window_inds], 'upr'],
          exposed = tag$exposed[window_inds]
        ) 
        
        ggplot(df, 
               aes(x = time, y = depth, ymin = depth_lwr, ymax = depth_upr,
                   alpha = factor(exposed))) + 
          # dive trace
          geom_pointrange() + 
          geom_line() + 
          # restrict trace to pre-exposure
          scale_alpha_manual(values = c('0'=1,'1'=0)) + 
          guides(alpha = 'none') + 
          # cee time
          geom_vline(xintercept = tag$exposure_time, lty = 3) + 
          # surface line
          geom_hline(yintercept = deep_dive_depth, lty = 3) + 
          # 1-step CTMC predictive distribution
          geom_point(
            data = data.frame(
              time = df$time[which(df$exposed == 1)[1]],
              depth = template_bins$center,
              depth_lwr = template_bins$center - template_bins$halfwidth,
              depth_upr = template_bins$center + template_bins$halfwidth,
              exposed = 0,
              p = Pmat[
                tag$depth.bin[window_inds][tail(which(df$exposed == 0),1)],
              ]
            ),
            mapping = aes(fill = p, size = p),
            shape = 22
          ) + 
          guides(size = 'none') + 
          scale_fill_distiller('Prob.', palette = 'PuBu', direction = 1,
                                breaks = seq(from = 0, to = .5, by = .1)) + 
          scale_size_continuous(range = c(3,7)) + 
          # formatting
          scale_y_reverse('Depth (m)', limits = rev(c(0, 1500))) + 
          theme_few() + 
          theme(text = element_text(size = 24))
      }
      
      pre_len = 1
      
      pl93_fast_descent = plot_tag_segment(
        tag = raw_sattags[[zc93_ind]], pre_len = pre_len, post_len = 1/6, 
        mvmt_type = 'fast_descent'
      ) + ggtitle('Fast descent') + xlab('Time (UTC)')
      
      pl93_slow_descent = plot_tag_segment(
        tag = raw_sattags[[zc93_ind]], pre_len = pre_len, post_len = 1/6, 
        mvmt_type = 'slow_descent'
      ) + ggtitle('Slow descent') + xlab('')
      
      pl93_fast_ascent = plot_tag_segment(
        tag = raw_sattags[[zc93_ind]], pre_len = pre_len, post_len = 1/6, 
        mvmt_type = 'fast_ascent'
      ) + ggtitle('Fast ascent') + xlab('Time (UTC)')
      
      pl93_slow_ascent = plot_tag_segment(
        tag = raw_sattags[[zc93_ind]], pre_len = pre_len, post_len = 1/6, 
        mvmt_type = 'slow_ascent'
      ) + ggtitle('Slow ascent') + xlab('')
      
      
      # save file
      f = file.path('output', 'figures')
      dir.create(f, showWarnings = FALSE, recursive = TRUE)
      ggsave(
        ggarrange(pl93_slow_descent, pl93_slow_ascent,
                  pl93_fast_descent, pl93_fast_ascent,
                  nrow = 2, ncol = 2, 
                  common.legend = TRUE, legend = 'right'), 
        filename = file.path(f, paste(tar_name(), '.png', sep = '')),
        dpi = 'print', width = 16, height = 8
      )
      
    }),
  
  tar_target(
    name = stage_dtmc_illustration,
    command = {
      
      fit_marginalized_model = list(
        samples = file.path('output', 'mcmc', 'fit_marginalized_model'),
        package = file.path('output', 'mcmc', 'fit_marginalized_model', 
                            'nim_pkg.rds')
      )
      
      #
      # load, label, and merge posterior samples
      #
      
      mvSample_files = dir(
        path = fit_marginalized_model$samples, 
        pattern = 'mvSamples_[0-9]+', 
        full.names = TRUE
      )
      
      mvSample2_files = dir(
        path = fit_marginalized_model$samples, 
        pattern = 'mvSamples2_[0-9]+', 
        full.names = TRUE
      )
      
      samples = do.call(rbind, lapply(mvSample_files, readRDS))
      samples2 = do.call(rbind, lapply(mvSample2_files, readRDS))
      
      colnames(samples) = readRDS(dir(
        path = fit_marginalized_model$samples, 
        pattern = 'mvSamples_colnames',
        full.names = TRUE
      ))
      
      colnames(samples2) = readRDS(dir(
        path = fit_marginalized_model$samples, 
        pattern = 'mvSamples2_colnames',
        full.names = TRUE
      ))
      
      samples = cbind(samples, samples2)
      rm(samples2)
      
      nim_pkg = readRDS(fit_marginalized_model$package)
      
      # set burn-in
      burn = 1:(nrow(samples)*.1)
      
      # get top-level names and groupings of variables being sampled
      sampling_targets = colnames(samples)
      sampling_target_groups = gsub(
        pattern = '\\[.*\\]',
        replacement = '',
        x = sampling_targets
      )
      sampling_groups = unique(sampling_target_groups)
      
      #
      # re-build model so that we can use it to generate transition matrices
      #
      
      # need to re-link compiled functions when running target in parallel
      source(file.path('R', 'util', 'statespace_tools', 'statespace_tools.R'))
      
      # uncompiled model
      mod = nimbleModel(
        code = modelCode, constants = nim_pkg$consts, data = nim_pkg$data,
        inits = nim_pkg$inits, name = tar_name(), calculate = FALSE
      )
      
      # compile model
      cmod = compileNimble(mod)
      
      # verify model has a finite likelihood
      cmod$calculate()
      
      #
      # plot an example progression of stages
      #
      
      tagName = 'ZcTag093'
      
      # determine subject id number for this tag
      subject_id = which(
        nim_pkg$consts$subject_id_labels == tagName
      )
      
      # get info. for the subject's final segment that was analyzed
      final_seg = nim_pkg$consts$segments[
        max(which(nim_pkg$consts$segments[,'subject_id'] == subject_id)),
      ]
      
      # transfer posterior mean of model parameters to model object
      for(tgt_group in sampling_groups) {
        cmod[[tgt_group]] = colMeans(samples[
          -burn, sampling_targets[sampling_target_groups == tgt_group]
        ])
      }
      
      # update model components
      cmod$calculate()
      
      # indices of the final segment
      seg_inds = final_seg['start_ind']:final_seg['end_ind']
      
      # ad-hoc estimate of latent state at last handful of observation times
      xf = do.call(c, lapply(12:1, function(lag) {
        pred_inds = seg_inds[1:(length(seg_inds)-lag+1)]
        p = finalPred2LayerCompressed(
          obs_lik_dict = cmod$depth_tx_mat,
          obs = cmod$depth_bins[pred_inds] - 1,
          txmat_dict = cmod$stage_tx_mat[final_seg['subject_id'], , , ],
          txmat_seq = cmod$covariateId[pred_inds] - 1,
          x0 = cmod$init_stages,
          num_obs_states = nim_pkg$consts$n_bins
        )
        which.max(p)
      }))
      
      # level order
      o = match(
        x = rev(c('fast_ascent', 'slow_ascent','fast_descent', 'slow_descent','slow_free')), 
        table = rownames(nim_pkg$consts$stage_defs)
      )
      
      df = data.frame(
        state = rownames(nim_pkg$consts$stage_defs)[xf],
        depth.bin = tail(nim_pkg$data$depth_bins[seg_inds], length(xf)),
        time = as.POSIXct(tail(data_pkg$data$times[seg_inds], length(xf)),
                          origin = '1970-01-01 00:00.00 UTC', tz = 'UTC')
      ) %>% mutate(
        state = factor(
          x = state, 
          labels = gsub('_', ' ', 
                        str_to_title(rownames(nim_pkg$consts$stage_defs))
                   )[o],
          levels = rownames(nim_pkg$consts$stage_defs)[o]
        ),
        depth = template_bins$center[depth.bin],
        lwr = depth - template_bins$halfwidth[depth.bin],
        upr = depth + template_bins$halfwidth[depth.bin]
      )
      
      df = rbind(df, data.frame(
        state = NA, time = tail(df$time,1) + 300, 
        depth = NA, lwr = NA, upr = NA, depth.bin = NA
      ))
      
      zc93_ind = which(sapply(raw_sattags, function(res) res$tag) == 'ZcTag093')
      tag = raw_sattags[[zc93_ind]]
      
      pl_mvmt = ggplot(df %>% filter(!is.na(state)), 
                       aes(x = time, y = state, group = 1)) + 
        # sequence of movement types
        geom_point() + 
        geom_line(lty = 3, col = 'grey60') + 
        # cee time
        geom_vline(xintercept = tag$exposure_time, lty = 3) + 
        # formatting
        scale_y_discrete('Mvmt. type', drop = FALSE) + 
        scale_x_datetime('Time (UTC)', limits = range(df$time)) +
        guides(x = 'none') + 
        theme_few() + 
        theme(axis.title.x = element_blank()) + 
        theme(text = element_text(size = 20))
      
      pl_depth = ggplot(df, aes(x = time, y = depth, ymin = lwr, ymax = upr)) + 
        # depth sequence
        geom_pointrange() + 
        geom_line() + 
        # cee time
        geom_vline(xintercept = tag$exposure_time, lty = 3) + 
        # surface 
        geom_hline(yintercept = template_bins$center[1], lty = 3) + 
        # formatting
        scale_x_datetime('Time (UTC)', limits = range(df$time)) +
        scale_y_reverse('Depth (m)', limits = rev(c(0, 1500))) + 
        theme_few() + 
        theme(text = element_text(size = 20))
        
      # save file
      f = file.path('output', 'figures')
      dir.create(f, showWarnings = FALSE, recursive = TRUE)
      ggsave(
        egg::ggarrange(pl_mvmt, pl_depth, nrow = 2), 
        filename = file.path(f, paste(tar_name(), '.png', sep = '')),
        dpi = 'print', width = 8, height = 4
      )
      
    }),
  
  sattag_illustration_script,
  
  tar_target(name = multiple_start_reps, command = 0:50),
  
  parameter_interpretation_plot_script,
  
  parameter_interpretation_variability_plot_script,
  
  parameter_table_script,
  
  parameter_interpretation_sensitivity_plot_script,
  
  random_effect_plot_script,
  
  random_starts_plot_script,
  
  random_starts_convergence_script
  
)
